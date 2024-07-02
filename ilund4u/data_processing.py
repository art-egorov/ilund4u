"""
This module contains all key data processing classes and methods of the tool.
"""
import pandas as pd
import numpy as np
import progress.bar
import collections
import subprocess
import typing
import shutil
import json
import time
import sys
import os

import leidenalg
import igraph

import Bio.Align.AlignInfo
import Bio.Data.CodonTable
import Bio.SeqRecord
import Bio.GenBank
import Bio.AlignIO
import Bio.Align
import Bio.SeqIO
import Bio.Seq
import BCBio.GFF

import msa4u

import ilund4u.manager
import ilund4u.methods
import ilund4u.drawing


class Proteome:
    """Proteome object represents a particular annotated proteome and its properties

    Attributes:
        proteome_id (str): Proteome identifier.
        gff_file (str): Path to the corresponding gff file.
        circular (int): [1 or 0] int value whether locus is circular or not. If genome is circular then first and last
            genes are considered as neighbours.
        cdss (pd.Series): List of annotated proteins.
        islands (pd.Series): Series of annotated islands.

    """

    def __init__(self, proteome_id: str, gff_file: str, circular: int = 1, cdss: typing.Union[None, pd.Series] = None,
                 islands: typing.Union[None, pd.Series] = None):
        """Proteome class constructor.

        Arguments:
            proteome_id (str): Proteome identifier.
            gff_file (str): Path to the corresponding gff file.
            circular (int): [1 or 0] int value whether locus is circular or not.
            cdss (pd.Series): List of annotated proteins.
            islands (pd.Series): Series of annotated islands.

        """
        self.proteome_id = proteome_id
        self.gff_file = gff_file
        self.cdss = cdss
        self.circular = circular
        self.islands = islands

    def get_proteome_db_row(self) -> dict:
        """Database building method for saving object's attributes.

        Returns:
            dict: modified object's attributes.

        """
        attributes_to_ignore = ["prms", "cdss", "islands"]
        attributes = {k: v for k, v in self.__dict__.items() if k not in attributes_to_ignore}
        attributes["gff_file"] = os.path.basename(attributes["gff_file"])
        attributes["cdss"] = self.cdss.apply(lambda cds: cds.cds_id).to_list()
        attributes["islands"] = self.islands.apply(lambda island: island.island_id).to_list()
        return attributes

    def annotate_variable_islands(self, prms: ilund4u.manager.Parameters) -> None:
        """Annotate proteome variable islands defined as a region with a set of non-conserved proteins.

        Arguments:
            prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        Returns:
            None

        """
        try:
            proteome_size = len(self.cdss.index)
            proteome_cds_classes = self.cdss.apply(lambda cds: cds.g_class).to_list()
            proteome_not_conserved_indexes = [ind for ind in range(proteome_size) if
                                              proteome_cds_classes[ind] != "conserved"]
            proteome_conserved_indexes = [ind for ind in range(proteome_size) if
                                          proteome_cds_classes[ind] == "conserved"]
            n_conserved_proteins = len(proteome_conserved_indexes)
            var_regions, cur_region = [], []
            self.islands = pd.Series()
            for pnci in proteome_not_conserved_indexes:
                if cur_region:
                    if pnci == cur_region[-1] + 1:
                        cur_region.append(pnci)
                    else:
                        var_regions.append(cur_region)
                        cur_region = [pnci]
                else:
                    cur_region.append(pnci)
                if pnci == proteome_not_conserved_indexes[-1] and cur_region:
                    var_regions.append(cur_region)
            if not var_regions:
                return None
            if self.circular:
                if var_regions[0][0] == 0 and var_regions[-1][-1] == proteome_size - 1:
                    last_region = var_regions.pop()
                    var_regions[0] = last_region + var_regions[0]
            proteome_islands_l = []
            for region in var_regions:
                var_indexes = [ind for ind in region if proteome_cds_classes[ind] == "variable"]
                if not var_indexes:
                    continue
                left_border, right_border = region[0], region[-1]
                left_cons_neighbours, right_cons_neighbours = [], []
                for dist in range(1, min(prms.args["neighbours_max_distance"], n_conserved_proteins // 2) + 1):
                    if self.circular:
                        left_index = (left_border - dist) % proteome_size
                        right_index = (right_border + dist) % proteome_size
                    else:
                        left_index = left_border - dist
                        right_index = right_border + dist
                    if left_index in proteome_conserved_indexes and \
                            len(left_cons_neighbours) < prms.args["neighbours_one_side_max_size"]:
                        left_cons_neighbours.append(left_index)
                    if right_index in proteome_conserved_indexes and \
                            len(right_cons_neighbours) < prms.args["neighbours_one_side_max_size"]:
                        right_cons_neighbours.append(right_index)
                left_cons_neighbours = left_cons_neighbours[::-1]
                cons_neighbours = left_cons_neighbours + right_cons_neighbours
                if len(set(cons_neighbours)) < prms.args["neighbours_min_size"]:
                    continue
                island_size = len(var_indexes)
                island_center = var_indexes[int(island_size / 2)]
                island_id = f"{self.proteome_id}:{island_center}"
                island = Island(island_id=island_id, proteome=self.proteome_id, center=island_center,
                                indexes=region, var_indexes=var_indexes, left_cons_neighbours=left_cons_neighbours,
                                right_cons_neighbours=right_cons_neighbours)
                proteome_islands_l.append(island)
            self.islands = pd.Series(proteome_islands_l, index=[pi.island_id for pi in proteome_islands_l])
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError(f"Unable to annotate variable islands for {self.proteome_id}") from error


class CDS:
    """CDS object represents an annotated protein.

    Attributes:
        cds_id (str): CDS identifier.
        proteome_id (str): Proteome identifier where CDS is encoded.
        start (int): 1-based start genomic coordinate.
        end (int): 1-based end genomic coordinates.
        length (int): length of the CDS.
        strand (int): Genomic strand (1: plus strand, -1: minus strand).
        name (str): Name of the feature which will be used as a label.
        group (str): CDS group that represents a set of homologous proteins.
        g_class (str): Class of CDS group (variable, intermediate, conserved).
        hmmscan_results (dict): Results of pyhmmer hmmscan annotation.

    """

    def __init__(self, cds_id: str, proteome_id: str, start: int, end: int, strand: int, name: str,
                 group: typing.Union[None, str] = None, g_class: typing.Union[None, str] = None,
                 hmmscan_results: typing.Union[None, dict] = None):
        """CDS class constructor.

        Arguments:
            cds_id (str): CDS identifier.
            proteome_id (str): Proteome identifier where CDS is encoded.
            start (int): 1-based start genomic coordinate.
            end (int): 1-based end genomic coordinates.
            strand (int): Genomic strand (1: plus strand, -1: minus strand).
            name (str): Name of the feature which will be used as a label.
            group (str): CDS group that represents a set of homologous proteins.
            g_class (str): Class of CDS group (variable, intermediate, conserved).
            hmmscan_results (dict): Results of pyhmmer hmmscan annotation.

        """
        self.cds_id = cds_id
        self.proteome_id = proteome_id
        self.start = start
        self.end = end
        self.length = int((end - start + 1) / 3)
        self.strand = strand
        self.name = name
        self.group = group
        self.g_class = g_class
        self.hmmscan_results = hmmscan_results

    def get_cds_db_row(self) -> dict:
        """Database building method for saving object's attributes.

        Returns:
            dict: object's attributes.

        """
        attributes_to_ignore = ["length"]
        attributes = {k: v for k, v in self.__dict__.items() if k not in attributes_to_ignore}
        return attributes


class Island:
    """Island object represents an annotated island defined as a region with a set of non-conserved proteins.

    Attributes:
        island_id (str): Island identifier.
        proteome (str): Proteome identifier where island is annotated.
        center (int): CDS index of the island center.
        indexes (list): CDS indexes of the island.
        size (int): Length of the island (number of CDSs).
        var_indexes (list): Indexes of CDS which g_class is "variable".
        left_cons_neighbours (list): Indexes of conserved neighbours on the left.
        right_cons_neighbours (list): Indexes of conserved neighbours on the right.
        databases_hits_stat (dict): Statistics from hmmscan annotation.

    """

    def __init__(self, island_id: str, proteome: str, center: int, indexes: list, var_indexes: list,
                 left_cons_neighbours: list, right_cons_neighbours: list,
                 databases_hits_stat: typing.Union[None, dict] = None):
        """Island class constructor.

        Arguments:
            island_id (str): Island identifier.
            proteome (str): Proteome identifier where island is annotated.
            center (int): CDS index of the island center.
            indexes (list): CDS indexes of the island.
            var_indexes (list): Indexes of CDS which g_class is "variable".
            left_cons_neighbours (list): Indexes of conserved neighbours on the left.
            right_cons_neighbours (list): Indexes of conserved neighbours on the right.
            databases_hits_stat (dict): Statistics from hmmscan annotation.

        """
        self.island_id = island_id
        self.proteome = proteome
        self.center = center
        self.indexes = indexes
        self.size = len(indexes)
        self.var_indexes = var_indexes
        self.left_cons_neighbours = left_cons_neighbours
        self.right_cons_neighbours = right_cons_neighbours
        if databases_hits_stat is None:
            databases_hits_stat = collections.defaultdict(lambda: collections.defaultdict(dict))
        self.databases_hits_stat = databases_hits_stat

    def get_island_db_row(self) -> dict:
        """Database building method for saving object's attributes.

        Returns:
            dict: object's attributes.

        """
        attributes_to_ignore = ["size"]
        attributes = {k: v for k, v in self.__dict__.items() if k not in attributes_to_ignore}
        return attributes

    def get_cons_neighbours_groups(self, cdss: pd.Series) -> list:
        """Get homology group attribute of conserved neighbour CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: Groups of conserved neighbours.

        """
        all_cons_neighbours = self.left_cons_neighbours + self.right_cons_neighbours
        cons_neighbours_groups = cdss.iloc[all_cons_neighbours].apply(lambda cds: cds.group).to_list()
        return cons_neighbours_groups

    def get_locus_groups(self, cdss: pd.Series) -> list:
        """Get homology group attribute of island locus CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: Groups of island CDSs.

        """
        locus_cdss_indexes = self.get_all_locus_indexes(cdss)
        locus_groups = cdss.iloc[locus_cdss_indexes].apply(lambda cds: cds.group).to_list()
        return locus_groups

    def get_island_groups(self, cdss: pd.Series) -> list:
        """Get homology group attribute of island CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: Groups of island CDSs.

        """
        island_groups = cdss.iloc[self.indexes].apply(lambda cds: cds.group).to_list()
        return island_groups

    def get_flanking_groups(self, cdss: pd.Series) -> list:
        """Get homology group attribute of island flanking CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: Groups of flanking CDSs.

        """
        locus_cdss_indexes = self.get_all_locus_indexes(cdss)
        flanking_indexes = [ind for ind in locus_cdss_indexes if ind not in self.indexes]
        flanking_groups = cdss.iloc[flanking_indexes].apply(lambda cds: cds.group).to_list()
        return flanking_groups

    def get_locus_proteins(self, cdss: pd.Series) -> list:
        """Get ids of locus CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: ids of island locus CDSs.

        """
        locus_cdss_indexes = self.get_all_locus_indexes(cdss)
        locus_groups = cdss.iloc[locus_cdss_indexes].apply(lambda cds: cds.cds_id).to_list()
        return locus_groups

    def get_island_proteins(self, cdss: pd.Series) -> list:
        """Get ids of island CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: Ids of island CDSs.

        """
        island_proteins = cdss.iloc[self.indexes].apply(lambda cds: cds.cds_id).to_list()
        return island_proteins

    def get_all_locus_indexes(self, cdss: pd.Series) -> list:
        """Get indexes of island locus CDSs.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            list: Indexes of island locus CDSs.

        """
        if self.left_cons_neighbours:
            island_left_border_cds_ind = self.left_cons_neighbours[0]
        else:
            island_left_border_cds_ind = self.indexes[0]
        if self.right_cons_neighbours:
            island_right_border_cds_ind = self.right_cons_neighbours[-1]
        else:
            island_right_border_cds_ind = self.indexes[-1]
        if island_left_border_cds_ind < island_right_border_cds_ind:
            island_cdss_indexes = [i for i in range(island_left_border_cds_ind, island_right_border_cds_ind + 1)]
        else:
            island_cdss_indexes = [i for i in range(island_left_border_cds_ind, cdss.size)] + \
                                  [i for i in range(0, island_right_border_cds_ind + 1)]
        return island_cdss_indexes

    def calculate_database_hits_stat(self, cdss: pd.Series) -> None:
        """Update statistics of hmmscan search results.

        Arguments:
            cdss (pd.Series): Series of annotated proteome CDSs.

        Returns:
            None

        """
        all_locus_indexes = self.get_all_locus_indexes(cdss)
        for ind in all_locus_indexes:
            lind_cds = cdss.iat[ind]
            lind_cds_hmmscan_res = lind_cds.hmmscan_results
            if lind_cds_hmmscan_res:
                if ind not in self.indexes:
                    self.databases_hits_stat[lind_cds_hmmscan_res["db"]]["flanking"][lind_cds.group] = \
                        lind_cds_hmmscan_res["target"]
                else:
                    self.databases_hits_stat[lind_cds_hmmscan_res["db"]]["cargo"][lind_cds.group] = \
                        lind_cds_hmmscan_res["target"]
        return None


class Proteomes:
    """Proteomes object represents a set of annotated genomes.

    Attributes:
        proteomes (pd.Series): Series (list) of Proteome objects.
        annotation (pd.DataFrame): Annotation table with description and statistics of proteomes.
        seq_to_ind (dict): Dictionary with proteome id to proteome index pairs.
        communities (dict): Dictionary representing annotated communities of proteomes
            (key - community_id, value - list of proteome ids)
        communities_annot (pd.DataFrame): Annotation table of proteome communities.
        proteins_fasta_file (str): Path to a fasta file containing all protein sequences.
        prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

    """

    def __init__(self, parameters: ilund4u.manager.Parameters):
        """Proteomes class constructor.

        Arguments:
            parameters (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        """
        self.proteomes = pd.Series()
        self.annotation = None
        self.__col_to_ind = None
        self.seq_to_ind = None
        self.communities = dict()
        self.communities_annot = None
        self.proteins_fasta_file = os.path.join(parameters.args["output_dir"], "all_proteins.fa")
        self.prms = parameters

    def save_as_db(self, db_folder: str) -> None:
        """Save Proteomes to the iLnd4u database.

        Arguments:
            db_folder (str): Database folder path.

        Returns:
            None

        """
        try:
            attributes_to_ignore = ["proteomes", "annotation", "communities_annot", "prms"]
            attributes = {k: v for k, v in self.__dict__.items() if k not in attributes_to_ignore}
            attributes["proteins_fasta_file"] = os.path.basename(attributes["proteins_fasta_file"])
            with open(os.path.join(db_folder, "proteomes.attributes.json"), 'w') as json_file:
                json.dump(attributes, json_file)
            self.communities_annot.to_csv(os.path.join(db_folder, "proteomes.communities_annot.tsv"), sep="\t",
                                          index_label="id")
            self.annotation["protein_clusters"] = self.annotation["protein_clusters"].apply(lambda x: ";".join(x))
            self.annotation.to_csv(os.path.join(db_folder, "proteomes.annotations.tsv"), sep="\t",
                                   index_label="proteome_id")
            os.mkdir(os.path.join(db_folder, "gff"))
            proteome_db_ind, cdss_db_ind, islands_db_ind, cds_ids, repr_cds_ids = [], [], [], [], set()
            for community, proteomes in self.communities.items():
                for proteome_id in proteomes:
                    proteome = self.proteomes.at[proteome_id]
                    proteome_db_ind.append(proteome.get_proteome_db_row())
                    os.system(f"cp '{proteome.gff_file}' {os.path.join(db_folder, 'gff')}/")
                    for cds in proteome.cdss.to_list():
                        cds_ids.append(cds.cds_id)
                        cdss_db_ind.append(cds.get_cds_db_row())
                        repr_cds_ids.add(cds.group)
                    for island in proteome.islands.to_list():
                        islands_db_ind.append(island.get_island_db_row())

            with open(os.path.join(db_folder, "proteome.ind.attributes.json"), "w") as json_file:
                json.dump(proteome_db_ind, json_file)
            with open(os.path.join(db_folder, "cds.ind.attributes.json"), "w") as json_file:
                json.dump(cdss_db_ind, json_file)
            with open(os.path.join(db_folder, "island.ind.attributes.json"), "w") as json_file:
                json.dump(islands_db_ind, json_file)

            initial_fasta_file = Bio.SeqIO.index(self.proteins_fasta_file, "fasta")
            with open(os.path.join(db_folder, attributes["proteins_fasta_file"]), "wb") as out_handle:
                for acc in cds_ids:
                    out_handle.write(initial_fasta_file.get_raw(acc))

            with open(os.path.join(db_folder, "representative_seqs.fa"), "wb") as out_handle:
                for acc in repr_cds_ids:
                    out_handle.write(initial_fasta_file.get_raw(acc))

            mmseqs_db_folder = os.path.join(db_folder, "mmseqs_db")
            if os.path.exists(mmseqs_db_folder):
                shutil.rmtree(mmseqs_db_folder)
            os.mkdir(mmseqs_db_folder)
            subprocess.run([self.prms.args["mmseqs_binary"], "createdb",
                            os.path.join(db_folder, attributes["proteins_fasta_file"]),
                            os.path.join(mmseqs_db_folder, "all_proteins")],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError(f"Unable to write Proteomes to the database.") from error

    @classmethod
    def db_init(cls, db_path: str, parameters: ilund4u.manager.Parameters):
        """Class method to load a Proteomes object from a database.

        Arguments:
            db_path (str): path to the database.
            parameters (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        Returns:
            cls: Proteomes object.

        """
        try:
            if parameters.args["verbose"]:
                print(f"○ Loading cds objects...", file=sys.stdout)
            with open(os.path.join(db_path, "cds.ind.attributes.json"), "r") as json_file:
                cds_ind_attributes = json.load(json_file)
            if parameters.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=len(cds_ind_attributes), suffix='%(index)d/%(max)d')
            cds_list = []
            for cds_dict in cds_ind_attributes:
                if parameters.args["verbose"]:
                    bar.next()
                cds_list.append(CDS(**cds_dict))
            if parameters.args["verbose"]:
                bar.finish()
            cdss = pd.Series(cds_list, index=[cds.cds_id for cds in cds_list])
            if parameters.args["verbose"]:
                print(f"○ Loading island objects...", file=sys.stdout)
            with open(os.path.join(db_path, "island.ind.attributes.json"), "r") as json_file:
                island_ind_attributes = json.load(json_file)
            if parameters.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=len(island_ind_attributes), suffix='%(index)d/%(max)d')
            island_list = []
            for proteome_dict in island_ind_attributes:
                if parameters.args["verbose"]:
                    bar.next()
                island_list.append(Island(**proteome_dict))
            if parameters.args["verbose"]:
                bar.finish()
            islands = pd.Series(island_list, index=[island.island_id for island in island_list])

            if parameters.args["verbose"]:
                print(f"○ Loading proteome objects...", file=sys.stdout)
            with open(os.path.join(db_path, "proteome.ind.attributes.json"), "r") as json_file:
                proteome_ind_attributes = json.load(json_file)
            if parameters.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=len(proteome_ind_attributes), suffix='%(index)d/%(max)d')
            proteome_list = []
            for proteome_dict in proteome_ind_attributes:
                if parameters.args["verbose"]:
                    bar.next()
                proteome_dict["gff_file"] = os.path.join(db_path, "gff", proteome_dict["gff_file"])
                proteome_dict["cdss"] = cdss.loc[proteome_dict["cdss"]]
                proteome_dict["islands"] = islands.loc[proteome_dict["islands"]]
                proteome_list.append(Proteome(**proteome_dict))
            if parameters.args["verbose"]:
                bar.finish()
            proteomes = pd.Series(proteome_list, index=[proteome.proteome_id for proteome in proteome_list])

            cls_obj = cls(parameters)
            cls_obj.proteomes = proteomes
            with open(os.path.join(db_path, "proteomes.attributes.json"), "r") as json_file:
                proteomes_attributes = json.load(json_file)
            cls_obj.communities = {int(k): v for k, v in proteomes_attributes["communities"].items()}

            cls_obj.proteins_fasta_file = os.path.join(db_path, proteomes_attributes["proteins_fasta_file"])

            cls_obj.annotation = pd.read_table(os.path.join(db_path, "proteomes.annotations.tsv"),
                                               sep="\t").set_index("proteome_id")
            cls_obj.__col_to_ind = {col: idx for idx, col in enumerate(cls_obj.annotation.columns)}
            cls_obj.communities_annot = pd.read_table(os.path.join(db_path, "proteomes.communities_annot.tsv"),
                                                      sep="\t").set_index("id")
            cls_obj.seq_to_ind = {sid: idx for idx, sid in enumerate(cls_obj.annotation.index)}

            return cls_obj
        except Exception as error:
            raise ilund4u.manager.ilund4uError(f"Unable to read Proteomes from the database.") from error

    def load_sequences_from_extended_gff(self, input_f: typing.Union[str, list], genome_annotation=None) -> None:
        """Load proteomes from gff files.

        Arguments:
            input_f (str | list): List of file paths or path to a folder with gff files.
            genome_annotation (path): Path to a table with annotation of genome circularity.
                Format: two columns with names: id, circular; tab-separated, 1,0 values.

        Returns:
            None

        """
        try:
            if isinstance(input_f, str):
                input_folder = input_f
                if not os.path.exists(input_folder):
                    raise ilund4u.manager.ilund4uError(f"Folder {input_folder} does not exist.")
                gff_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if
                             os.path.splitext(f)[-1].lower() == ".gff"]
            elif isinstance(input_f, list):
                gff_files = input_f
            else:
                raise ilund4u.manager.ilund4uError(f"The input for the GFF parsing function must be either a folder or "
                                                   f"a list of files.")
            if not gff_files:
                raise ilund4u.manager.ilund4uError(f"Folder {input_f} does not contain .gff(.GFF) files.")
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            else:
                if os.path.exists(self.proteins_fasta_file):
                    os.remove(self.proteins_fasta_file)
            genome_circularity_dict = dict()
            if genome_annotation:
                try:
                    genome_annotation_table = pd.read_table(genome_annotation, sep="\t").set_index("id")
                    genome_circularity_dict = genome_annotation_table["circular"].to_dict()
                except:
                    print("○ Warning: unable to read genome annotation table. Check the format.")
            num_of_gff_files = len(gff_files)
            if self.prms.args["verbose"]:
                print(f"○ Reading gff file{'s' if len(gff_files) > 1 else ''}...", file=sys.stdout)
            if len(gff_files) > 1 and self.prms.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=num_of_gff_files, suffix='%(index)d/%(max)d')
            proteome_list, annotation_rows = [], []
            gff_records_batch = []
            for gff_file_index, gff_file_path in enumerate(gff_files):
                try:
                    if len(gff_files) > 1 and self.prms.args["verbose"]:
                        bar.next()
                    gff_records = list(BCBio.GFF.parse(gff_file_path, limit_info=dict(gff_type=["CDS"])))
                    if len(gff_records) != 1:
                        print(f"\n○ Warning: gff file {gff_file_path} contains information for more than 1 "
                              f"sequence. File will be skipped.")
                        continue
                    gff_record = gff_records[0]
                    try:
                        record_locus_sequence = gff_record.seq
                    except Bio.Seq.UndefinedSequenceError as error:
                        raise ilund4u.manager.ilund4uError(f"gff file doesn't contain corresponding "
                                                           f"sequences.") from error
                    features_ids = [i.id for i in gff_record.features]
                    if len(features_ids) != len(set(features_ids)):
                        raise ilund4u.manager.ilund4uError(f"Gff file {gff_file_path} contains duplicated feature "
                                                           f"ids while only unique are allowed.")
                    if len(features_ids) > self.prms.args["min_proteome_size"]:
                        if gff_record.id in genome_circularity_dict.keys():
                            circular = int(genome_circularity_dict[gff_record.id])
                        else:
                            circular = int(self.prms.args["circular_genomes"])
                        record_proteome = Proteome(proteome_id=gff_record.id, gff_file=gff_file_path, cdss=pd.Series(),
                                                   circular=circular)
                        record_cdss = []
                        gff_records = []
                        for gff_feature in gff_record.features:
                            cds_id = gff_feature.id.replace(";", ",")
                            if gff_record.id not in cds_id:
                                cds_id = f"{gff_record.id}-{cds_id}"  # Attention
                            transl_table = self.prms.args["default_transl_table"]
                            if "transl_table" in gff_feature.qualifiers.keys():
                                transl_table = int(gff_feature.qualifiers["transl_table"][0])
                            name = ""
                            if self.prms.args["gff_CDS_name_source"] in gff_feature.qualifiers:
                                name = gff_feature.qualifiers[self.prms.args["gff_CDS_name_source"]][0]

                            sequence = gff_feature.translate(record_locus_sequence, table=transl_table, cds=False)[:-1]
                            gff_records_batch.append(Bio.SeqRecord.SeqRecord(seq=sequence, id=cds_id, description=""))
                            cds = CDS(cds_id=cds_id, proteome_id=gff_record.id,
                                      start=int(gff_feature.location.start) + 1, end=int(gff_feature.location.end),
                                      strand=gff_feature.location.strand, name=name)
                            record_cdss.append(cds)
                        record_proteome.cdss = pd.Series(record_cdss, index=[cds.cds_id for cds in record_cdss])
                        proteome_list.append(record_proteome)
                        annotation_rows.append(dict(id=gff_record.id, length=len(gff_record.seq),
                                                    proteome_size=len(features_ids),
                                                    proteome_size_unique="", protein_clusters=""))
                    if gff_file_index % 1000 == 0 or gff_file_index == len(gff_files) - 1:
                        with open(self.proteins_fasta_file, "a") as handle:
                            Bio.SeqIO.write(gff_records_batch, handle, "fasta")
                        gff_records_batch = []
                except:
                    print(f"○ Warning: gff file {gff_file_path} was not read properly")
            if len(gff_files) > 1 and self.prms.args["verbose"]:
                bar.finish()
            self.proteomes = pd.Series(proteome_list, index=[pr.proteome_id for pr in proteome_list])
            self.annotation = pd.DataFrame(annotation_rows).set_index("id")
            self.__col_to_ind = {col: idx for idx, col in enumerate(self.annotation.columns)}
            self.annotation = self.annotation.sort_values(by="proteome_size")
            self.proteomes = self.proteomes.loc[self.annotation.index]
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to load proteomes from gff files.") from error

    def mmseqs_cluster(self) -> dict:
        """Cluster all proteins using mmseqs in order to define groups of homologues.

        Returns:
            dict: protein id to cluster id dictionary.

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Running mmseqs for protein clustering...", file=sys.stdout)
            mmseqs_input = self.proteins_fasta_file
            mmseqs_output_folder = os.path.join(self.prms.args["output_dir"], "mmseqs")
            if os.path.exists(mmseqs_output_folder):
                shutil.rmtree(mmseqs_output_folder)
            os.mkdir(mmseqs_output_folder)
            mmseqs_output_folder_db = os.path.join(mmseqs_output_folder, "DB")
            os.mkdir(mmseqs_output_folder_db)
            mmseqs_stdout = open(os.path.join(mmseqs_output_folder, "mmseqs_stdout.txt"), "w")
            mmseqs_stderr = open(os.path.join(mmseqs_output_folder, "mmseqs_stderr.txt"), "w")
            subprocess.run([self.prms.args["mmseqs_binary"], "createdb", mmseqs_input,
                            os.path.join(mmseqs_output_folder_db, "sequencesDB")], stdout=mmseqs_stdout,
                           stderr=mmseqs_stderr)
            subprocess.run([self.prms.args["mmseqs_binary"], "cluster",
                            os.path.join(mmseqs_output_folder_db, "sequencesDB"),
                            os.path.join(mmseqs_output_folder_db, "clusterDB"),
                            os.path.join(mmseqs_output_folder_db, "tmp"),
                            "--cluster-mode", str(self.prms.args["mmseqs_cluster_mode"]),
                            "--cov-mode", str(self.prms.args["mmseqs_cov_mode"]),
                            "--min-seq-id", str(self.prms.args["mmseqs_min_seq_id"]),
                            "-c", str(self.prms.args["mmseqs_c"]),
                            "-s", str(self.prms.args["mmseqs_s"])], stdout=mmseqs_stdout,
                           stderr=mmseqs_stderr)  # threads!
            subprocess.run([self.prms.args["mmseqs_binary"], "createtsv",
                            os.path.join(mmseqs_output_folder_db, "sequencesDB"),
                            os.path.join(mmseqs_output_folder_db, "sequencesDB"),
                            os.path.join(mmseqs_output_folder_db, "clusterDB"),
                            os.path.join(mmseqs_output_folder, "mmseqs_clustering.tsv")],
                           stdout=mmseqs_stdout, stderr=mmseqs_stderr)
            mmseqs_clustering_results = pd.read_table(os.path.join(mmseqs_output_folder, "mmseqs_clustering.tsv"),
                                                      sep="\t", header=None, names=["cluster", "protein_id"])
            mmseqs_clustering_results = mmseqs_clustering_results.set_index("protein_id")["cluster"].to_dict()
            num_of_unique_clusters = len(set(mmseqs_clustering_results.values()))
            num_of_proteins = len(mmseqs_clustering_results.keys())
            if self.prms.args["verbose"]:
                print(f"  ⦿ {num_of_unique_clusters} clusters for {num_of_proteins} proteins were found with mmseqs\n"
                      f"  ⦿ mmseqs clustering results were saved to "
                      f"{os.path.join(mmseqs_output_folder, 'mmseqs_clustering.tsv')}", file=sys.stdout)
            return mmseqs_clustering_results
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to run mmseqs clustering.") from error

    def process_mmseqs_results(self, mmseqs_results: dict) -> dict:
        """Process results of mmseqs clustering run.

        Arguments:
            mmseqs_results (dict): results of mmseqs_cluster function.

        Returns:
            dict: dictionary with protein cluster id to list of protein ids items.

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Processing mmseqs results ...", file=sys.stdout)
            sequences_to_drop, drop_reason = [], []
            current_p_length, cpl_added_proteomes = None, None
            cluster_to_sequences = collections.defaultdict(list)
            if self.prms.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=len(self.proteomes.index), suffix='%(index)d/%(max)d')
            for p_index, proteome in enumerate(self.proteomes.to_list()):
                if self.prms.args["verbose"]:
                    bar.next()
                seq_p_size = self.annotation.iat[p_index, self.__col_to_ind["proteome_size"]]
                if seq_p_size != current_p_length:
                    current_p_length = seq_p_size
                    cpl_added_proteomes = []
                seq_protein_clusters = []
                for cds in proteome.cdss.to_list():
                    cds.group = mmseqs_results[cds.cds_id]
                    seq_protein_clusters.append(cds.group)
                seq_protein_clusters_set = set(seq_protein_clusters)
                unique_p_size = len(seq_protein_clusters_set)
                if seq_protein_clusters_set in cpl_added_proteomes:
                    sequences_to_drop.append(proteome.proteome_id)
                    drop_reason.append(f"Duplicate")
                    continue
                if unique_p_size / seq_p_size < self.prms.args["proteome_uniqueness_cutoff"]:
                    sequences_to_drop.append(proteome.proteome_id)
                    drop_reason.append("Proteome uniqueness cutoff")
                    continue
                self.annotation.iat[p_index, self.__col_to_ind["proteome_size_unique"]] = unique_p_size
                self.annotation.iat[p_index, self.__col_to_ind["protein_clusters"]] = list(seq_protein_clusters_set)
                cpl_added_proteomes.append(seq_protein_clusters_set)
                for p_cluster in seq_protein_clusters_set:
                    cluster_to_sequences[p_cluster].append(proteome.proteome_id)
            dropped_sequences = pd.DataFrame(dict(sequence=sequences_to_drop, reason=drop_reason))
            dropped_sequences.to_csv(os.path.join(self.prms.args["output_dir"], "dropped_sequences.tsv"), sep="\t",
                                     index=False)
            if self.prms.args["verbose"]:
                bar.finish()
                print(f"  ⦿ {len(sequences_to_drop)} proteomes were excluded after proteome"
                      f" deduplication and filtering", file=sys.stdout)
            self.annotation = self.annotation.drop(sequences_to_drop)
            self.proteomes = self.proteomes.drop(sequences_to_drop)
            self.annotation = self.annotation.sort_values(by="proteome_size_unique")
            self.proteomes = self.proteomes.loc[self.annotation.index]
            self.annotation["index"] = list(range(len(self.proteomes.index)))
            self.seq_to_ind = {sid: idx for idx, sid in enumerate(self.annotation.index)}
            return cluster_to_sequences
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to process mmseqs output.") from error

    def build_proteome_network(self, cluster_to_sequences: dict) -> igraph.Graph:
        """Build proteome network where each proteome represented by node and weighted edges between nodes -
            fraction of shared homologues.

        Arguments:
            cluster_to_sequences (dict): cluster id to list of proteins dictionary
                (results of process_mmseqs_results() function)

        Returns:
            igraph.Graph: proteome network.

        """
        try:
            cluster_to_proteome_index = dict()
            for cluster, sequences in cluster_to_sequences.items():
                indexes = sorted([self.seq_to_ind[seq_id] for seq_id in sequences])
                cluster_to_proteome_index[cluster] = collections.deque(indexes)

            proteome_sizes = self.annotation[["proteome_size_unique", "index"]]
            first_index_for_size = proteome_sizes.groupby("proteome_size_unique").tail(1).copy()
            max_p_size = first_index_for_size["proteome_size_unique"].max()
            cut_off_mult = 1 / self.prms.args["proteome_similarity_cutoff"]
            first_index_for_size["cutoff"] = first_index_for_size["proteome_size_unique"].apply(
                lambda size: first_index_for_size[first_index_for_size["proteome_size_unique"] >=
                                                  min(size * cut_off_mult, max_p_size)]["index"].min() + 1)
            upper_index_cutoff = first_index_for_size.set_index("proteome_size_unique")["cutoff"]
            proteome_sizes = proteome_sizes.set_index("index")["proteome_size_unique"]

            if self.prms.args["verbose"]:
                print(f"○ Proteomes network construction...", file=sys.stdout)
                bar = progress.bar.FillingCirclesBar(" ", max=len(self.proteomes.index), suffix="%(index)d/%(max)d")
            stime = time.time()
            edges, weights = [], []
            for i in range(len(self.proteomes.index)):
                if self.prms.args["verbose"]:
                    bar.next()
                clusters_i = self.annotation.iat[i, self.__col_to_ind["protein_clusters"]]
                size_i = self.annotation.iat[i, self.__col_to_ind["proteome_size_unique"]]
                counts_i = collections.defaultdict(int)
                upper_i_cutoff = upper_index_cutoff.at[size_i]
                for cl in clusters_i:
                    js = cluster_to_proteome_index[cl]
                    for j in js.copy():
                        if i < j < upper_i_cutoff:
                            counts_i[j] += 1
                        elif j <= i:
                            js.popleft()
                        else:
                            break
                weights_i = pd.Series(counts_i)
                proteome_sizes_connected = proteome_sizes.iloc[weights_i.index]
                norm_factor_i = pd.Series(
                    0.5 * (size_i + proteome_sizes_connected) / (size_i * proteome_sizes_connected), \
                    index=weights_i.index)
                weights_i = weights_i.mul(norm_factor_i)
                weights_i = weights_i[weights_i >= self.prms.args["proteome_similarity_cutoff"]]
                for j, w in weights_i.items():
                    edges.append([i, j])
                    weights.append(round(w, 4))
            if self.prms.args["verbose"]:
                bar.finish()
            etime = time.time()
            if self.prms.args["verbose"]:
                print(f"  ⦿ Network building elapsed time: {round(etime - stime, 2)} sec")
            graph = igraph.Graph(len(self.proteomes.index), edges, directed=False)
            graph.vs["index"] = self.annotation["index"].to_list()
            graph.vs["sequence_id"] = self.annotation.index.to_list()
            graph.es["weight"] = weights
            graph.save(os.path.join(self.prms.args["output_dir"], "proteome_network.gml"))
            if self.prms.args["verbose"]:
                print(f"  ⦿ Proteomes network with {len(edges)} connections was built\n"
                      f"  ⦿ Network was saved as {os.path.join(self.prms.args['output_dir'], 'proteome_network.gml')}",
                      file=sys.stdout)
            return graph
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to built proteome network.") from error

    def find_proteome_communities(self, graph: igraph.Graph) -> None:
        """Find proteome communities using Leiden algorithm and update communities attribute.

        Arguments:
            graph (igraph.Graph): network of proteomes obtained by build_proteome_network function.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print("○ Proteome network partitioning using the Leiden algorithm...")
            partition_leiden = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition,
                                                        resolution_parameter=self.prms.args[
                                                            "leiden_resolution_parameter_p"],
                                                        weights="weight", n_iterations=-1)
            graph.vs["communities_Leiden"] = partition_leiden.membership
            if self.prms.args["verbose"]:
                print(f"  ⦿ {len(set(partition_leiden.membership))} proteome communities were found")
            communities_annot_rows = []
            sequences_to_drop = []
            for community_index, community in enumerate(partition_leiden):
                community_size = len(community)
                subgraph = graph.subgraph(community)
                proteomes = subgraph.vs["sequence_id"]
                if community_size >= self.prms.args["min_proteome_community_size"]:
                    self.communities[community_index] = proteomes
                else:
                    sequences_to_drop += proteomes
                if community_size > 1:
                    subgraph_edges = subgraph.get_edgelist()
                    num_of_edges = len(subgraph_edges)
                    num_of_edges_fr = num_of_edges / (community_size * (community_size - 1) * 0.5)
                    weights = subgraph.es["weight"]
                    avg_weight = round(np.mean(weights), 3)
                    max_identity = max(weights)
                else:
                    num_of_edges, num_of_edges_fr, avg_weight, max_identity = "", "", "", ""
                communities_annot_rows.append([community_index, community_size, avg_weight, max_identity,
                                               num_of_edges_fr, num_of_edges, ";".join(proteomes)])
            communities_annot = pd.DataFrame(communities_annot_rows, columns=["id", "size", "avg_weight", "max_weight",
                                                                              "fr_edges", "n_edges", "proteomes"])
            communities_annot.to_csv(os.path.join(os.path.join(self.prms.args["output_dir"],
                                                               "proteome_communities.tsv")), sep="\t", index=False)
            communities_annot = communities_annot[communities_annot["size"] >=
                                                  self.prms.args["min_proteome_community_size"]]
            self.communities_annot = communities_annot.set_index("id")
            self.annotation = self.annotation.drop(sequences_to_drop)
            self.proteomes = self.proteomes.drop(sequences_to_drop)
            self.annotation["index"] = list(range(len(self.proteomes.index)))
            self.seq_to_ind = {sid: idx for idx, sid in enumerate(self.annotation.index)}
            if self.prms.args["verbose"]:
                print(f"  ⦿ {len(communities_annot.index)} proteomes communities with size >= "
                      f"{self.prms.args['min_proteome_community_size']} were taken for further analysis",
                      file=sys.stdout)
                print(f"  ⦿ {len(sequences_to_drop)} proteomes from smaller communities were excluded from the "
                      f"analysis", file=sys.stdout)
            if len(communities_annot.index) == 0:
                print("○ Termination since no proteome community was taken for further analysis")
                sys.exit()
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to find proteome communities.") from error

    def define_protein_classes(self) -> None:
        """Define protein classes (conserved, intermediate, variable) based on presence in a proteome community.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print("○ Defining protein classes within each community...")
            number_of_communities = len(self.communities_annot.index)
            protein_classes_trows = []
            if self.prms.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=number_of_communities, suffix='%(index)d/%(max)d')
            for com_id, com_pr_ids in self.communities.items():
                if self.prms.args["verbose"]:
                    bar.next()
                com_size = len(com_pr_ids)
                com_annotation = self.annotation.loc[com_pr_ids]
                com_protein_clusters = com_annotation["protein_clusters"].sum()
                com_protein_clusters_count = collections.Counter(com_protein_clusters)
                com_protein_classes = dict()
                for pc, counts in com_protein_clusters_count.items():
                    pc_fraction = counts / com_size
                    if pc_fraction < self.prms.args["variable_protein_cluster_cutoff"]:
                        pc_class = "variable"
                    elif pc_fraction > self.prms.args["conserved_protein_cluster_cutoff"]:
                        pc_class = "conserved"
                    else:
                        pc_class = "intermediate"
                    com_protein_classes[pc] = pc_class
                    protein_classes_trows.append(dict(community=com_id, community_size=com_size, protein_group=pc,
                                                      protein_group_class=pc_class, fraction=pc_fraction,
                                                      protein_group_counts=counts))
                com_proteomes = self.proteomes.loc[com_pr_ids]
                for com_proteome in com_proteomes:
                    for cds in com_proteome.cdss:
                        cds.g_class = com_protein_classes[cds.group]
            if self.prms.args["verbose"]:
                bar.finish()
            protein_classes_t = pd.DataFrame(protein_classes_trows)
            protein_classes_t.to_csv(os.path.join(self.prms.args["output_dir"], "protein_group_classes.tsv"), sep="\t",
                                     index=False)
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to define protein classes.") from error

    def annotate_variable_islands(self) -> None:
        """Annotate variable islands defined as a region with a set of non-conserved proteins.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Annotating variable islands within each proteome...", file=sys.stdout)
            total_number_of_variable_regions = 0
            for proteome_index, proteome in enumerate(self.proteomes):
                proteome.annotate_variable_islands(self.prms)
                total_number_of_variable_regions += proteome.islands.size
            if self.prms.args["verbose"]:
                print(f"  ⦿ {total_number_of_variable_regions} variable regions are annotated in "
                      f"{len(self.proteomes.index)} proteomes "
                      f"({round(total_number_of_variable_regions / len(self.proteomes.index), 3)} per proteome)",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to annotate variable islands.") from error

    def build_islands_network(self) -> igraph.Graph:
        """Build island network where each node - an island and weighted edges - fraction of shared
            conserved neighbours homologues.

        Returns:
            igraph.Graph: Island network.

        """
        try:
            if self.prms.args["verbose"]:
                print("○ Island network construction within each proteome community...")
            output_network_folder = os.path.join(self.prms.args["output_dir"], "island_networks")
            if os.path.exists(output_network_folder):
                shutil.rmtree(output_network_folder)
            os.mkdir(output_network_folder)
            number_of_communities = len(self.communities_annot.index)
            stime = time.time()
            networks = []
            if self.prms.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=number_of_communities, suffix="%(index)d/%(max)d")
            for com_id, com_pr_ids in self.communities.items():
                if self.prms.args["verbose"]:
                    bar.next()
                com_proteomes = self.proteomes.loc[com_pr_ids]
                com_island_n_sizes = pd.Series()
                com_neighbours = pd.Series()
                cluster_to_island = collections.defaultdict(collections.deque)

                islands_list = [island for proteome in com_proteomes.to_list() for island in proteome.islands.to_list()]
                island_id_to_index = {isl.island_id: ind for ind, isl in enumerate(islands_list)}
                for proteome in com_proteomes:
                    for island in proteome.islands.to_list():
                        island_id = island.island_id
                        island_index = island_id_to_index[island_id]
                        conserved_island_neighbours_groups = set(island.get_cons_neighbours_groups(proteome.cdss))
                        com_neighbours.at[island_index] = list(conserved_island_neighbours_groups)
                        com_island_n_sizes.at[island_index] = len(conserved_island_neighbours_groups)
                        for cing in conserved_island_neighbours_groups:
                            cluster_to_island[cing].append(island_index)
                edges, weights = [], []

                for i in range(len(com_island_n_sizes.index)):
                    neighbours_i = com_neighbours.iat[i]
                    size_i = com_island_n_sizes.iat[i]
                    counts_i = collections.defaultdict(int)
                    for ncl in neighbours_i:
                        js = cluster_to_island[ncl]
                        for j in js.copy():
                            if i < j:
                                counts_i[j] += 1
                            else:
                                js.popleft()
                    weights_i = pd.Series(counts_i)
                    connected_n_sizes = com_island_n_sizes.iloc[weights_i.index]
                    norm_factor_i = pd.Series(0.5 * (size_i + connected_n_sizes) / (size_i * connected_n_sizes), \
                                              index=weights_i.index)
                    weights_i = weights_i.mul(norm_factor_i)
                    weights_i = weights_i[weights_i >= self.prms.args["island_neighbours_similarity_cutoff"]]
                    for j, w in weights_i.items():
                        edges.append([i, j])
                        weights.append(round(w, 4))
                graph = igraph.Graph(len(com_island_n_sizes.index), edges, directed=False)
                graph.vs["index"] = com_island_n_sizes.index.to_list()
                graph.vs["island_id"] = [isl.island_id for isl in islands_list]
                graph.vs["island_size"] = [isl.size for isl in islands_list]
                graph.vs["proteome_id"] = [isl.proteome for isl in islands_list]
                graph.es["weight"] = weights
                graph.save(os.path.join(output_network_folder, f"{com_id}.gml"))
                networks.append(graph)
            if self.prms.args["verbose"]:
                bar.finish()
            etime = time.time()
            if self.prms.args["verbose"]:
                print(f"  ⦿ Island network building elapsed time: {round(etime - stime, 2)} sec")
            return networks
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to build island network.") from error

    def find_hotspots(self, networks: igraph.Graph) -> None:
        """Find hotspots in an island network using Leiden algorithm.

        Args:
            networks (igraph.Graph): Island network obtained by build_islands_network() function.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print("○ Searching for hotspots within each community...")
            number_of_communities = len(self.communities_annot.index)
            if self.prms.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=number_of_communities, suffix="%(index)d/%(max)d")
            hotspots_l, h_annotation_rows = [], []
            for com_id, com_pr_ids in self.communities.items():
                if self.prms.args["verbose"]:
                    bar.next()
                com_size = len(com_pr_ids)
                com_proteomes = self.proteomes.loc[com_pr_ids]
                com_island_network = networks[com_id]
                partition_leiden = leidenalg.find_partition(com_island_network, leidenalg.CPMVertexPartition,
                                                            resolution_parameter=self.prms.args[
                                                                "leiden_resolution_parameter_i"],
                                                            weights="weight", n_iterations=-1)
                com_island_network.vs["communities_Leiden"] = partition_leiden.membership
                islands_list = [island for proteome in com_proteomes for island in proteome.islands.to_list()]
                for icom_ind, i_com in enumerate(partition_leiden):
                    subgraph = com_island_network.subgraph(i_com)
                    proteomes = subgraph.vs["proteome_id"]
                    hotspot_uniq_size = len(set(proteomes))
                    hotspot_presence = hotspot_uniq_size / com_size
                    if hotspot_presence > self.prms.args["hotspot_presence_cutoff"]:
                        island_indexes, island_ids = subgraph.vs["index"], subgraph.vs["island_id"]
                        island_size = subgraph.vs["island_size"]
                        strength, degree = subgraph.strength(weights="weight"), subgraph.degree()
                        island_annotation = pd.DataFrame(dict(island=island_ids, island_index=island_indexes,
                                                              island_size=island_size,
                                                              proteome=proteomes, strength=strength,
                                                              degree=degree)).set_index("island")
                        if self.prms.args["deduplicate_proteomes_within_hotspot"]:  # To update usage wo it
                            island_annotation = island_annotation.sort_values(by="strength", ascending=False)
                            island_annotation = island_annotation.drop_duplicates(subset="proteome", keep="first")
                            island_annotation = island_annotation.sort_values(by="island_index")
                            nodes_to_remove = subgraph.vs.select(island_id_notin=island_annotation.index.to_list())
                            subgraph.delete_vertices(nodes_to_remove)
                        islands = [islands_list[ind] for ind in island_annotation["island_index"].to_list()]
                        unique_island_cds_groups = []
                        conserved_island_groups_count = collections.defaultdict(int)
                        for island in islands:
                            island_proteome_cdss = com_proteomes.at[island.proteome].cdss
                            island_cds_groups = island_proteome_cdss.iloc[island.indexes].apply(
                                lambda isl: isl.group).to_list()
                            ic_indexes = island.left_cons_neighbours + island.right_cons_neighbours
                            island_conserved_groups = island_proteome_cdss.iloc[ic_indexes].apply(
                                lambda isl: isl.group).to_list()
                            if set(island_cds_groups) not in unique_island_cds_groups:
                                unique_island_cds_groups.append(set(island_cds_groups))
                            for icg in set(island_conserved_groups):
                                conserved_island_groups_count[icg] += 1

                        signature_cutoff = int(self.prms.args["hotspot_signature_presence_cutoff"] * hotspot_uniq_size)
                        hotspot_conserved_signature = [g for g, c in conserved_island_groups_count.items() if
                                                       c >= signature_cutoff]

                        number_of_unique_islands = len(unique_island_cds_groups)
                        hotspot = Hotspot(hotspot_id=f"{com_id}-{icom_ind}", proteome_community=com_id,
                                          size=hotspot_uniq_size, islands=islands,
                                          conserved_signature=hotspot_conserved_signature,
                                          island_annotation=island_annotation)
                        h_annotation_row = dict(hotspot_id=f"{com_id}-{icom_ind}", size=hotspot_uniq_size,
                                                uniqueness=round(number_of_unique_islands / hotspot_uniq_size, 3),
                                                number_of_unique_islands=number_of_unique_islands,
                                                proteome_community=com_id)
                        h_annotation_rows.append(h_annotation_row)
                        hotspots_l.append(hotspot)
            if self.prms.args["verbose"]:
                bar.finish()
            h_annotation = pd.DataFrame(h_annotation_rows).set_index("hotspot_id")
            hotspots_s = pd.Series(hotspots_l, index=[hotspot.hotspot_id for hotspot in hotspots_l])
            hotspots_obj = Hotspots(hotspots_s, h_annotation, parameters=self.prms)
            num_of_hotspots = len(hotspots_l)
            if self.prms.args["verbose"]:
                print(f"  ⦿ {num_of_hotspots} hotspots were found in {number_of_communities} proteome communities"
                      f"  (Avg: {round(num_of_hotspots / number_of_communities, 3)} per community)",
                      file=sys.stdout)
            return hotspots_obj
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to find communities in the island network.") from error


class Hotspot:
    """Hostpot object represent a hotspot as a set of islands.

    Attributes:
        hotspot_id (str): Hotspot identifier.
        size (int): Number of islands.
        proteome_community (int): Identifier of proteome community where hotspot is annotated.
        islands (list): List of islands.
        conserved_signature (list): Conserved flanking proteins that are usually found in islands.
        island_annotation (pd.DataFrame): Annotation table of islands.

    """

    def __init__(self, hotspot_id: str, size: int, proteome_community: int, islands: list,
                 conserved_signature: list, island_annotation: pd.DataFrame):
        """Hotspot class constructor.

        Arguments:
            hotspot_id (str): Hotspot identifier.
            size (int): Number of islands.
            proteome_community (int): Identifier of proteome community where hotspot is annotated.
            islands (list): List of islands.
            conserved_signature (list): Conserved flanking proteins that are usually found in islands.
            island_annotation (pd.DataFrame): Annotation table of islands.

        """
        self.hotspot_id = hotspot_id
        self.size = size
        self.proteome_community = proteome_community
        self.islands = islands
        self.island_annotation = island_annotation
        self.conserved_signature = conserved_signature

    def get_hotspot_db_row(self) -> dict:
        """Database building method for saving object's attributes.

        Returns:
            dict: object's attributes.

        """
        attributes_to_ignore = ["island_annotation", "islands"]
        attributes = {k: v for k, v in self.__dict__.items() if k not in attributes_to_ignore}
        return attributes

    def calculate_database_hits_stats(self, proteomes: Proteomes,
                                      prms: ilund4u.manager.Parameters) -> collections.defaultdict:
        """Calculate statistics of pyhmmer annotation for island proteins.

        Arguments:
            proteomes (Proteomes): Proteomes object.
            prms (ilund4u.manager.Parameters): Parameters object.

        Returns:
            collections.defaultdict: hits to the databases.

        """
        hotspot_stat = collections.defaultdict(lambda: collections.defaultdict(dict))
        for island in self.islands:
            proteome = proteomes.proteomes.at[island.proteome]
            island.calculate_database_hits_stat(proteome.cdss)
            island_dbstat = island.databases_hits_stat
            db_names = prms.args["database_names"]
            for db_name in db_names:
                r_types = ["cargo", "flanking"]
                for r_type in r_types:
                    self.island_annotation.at[island.island_id, f"{db_name}_{r_type}"] = \
                        ",".join(island_dbstat[db_name][r_type].values())
                    hotspot_stat[db_name][r_type].update(island_dbstat[db_name][r_type])
        return hotspot_stat

    def get_hotspot_groups(self, proteomes: Proteomes) -> dict:
        """Get protein groups found on island cargo or as flanking genes

        Arguments:
            proteomes (Proteomes): Proteomes object.

        Returns:
            dict: cargo and flanking gene protein groups.

        """
        groups = dict(cargo=set(), flanking=set())
        for island in self.islands:
            proteome = proteomes.proteomes.at[island.proteome]
            groups["cargo"].update(set(island.get_island_groups(proteome.cdss)))
            groups["flanking"].update(set(island.get_flanking_groups(proteome.cdss)))
        return groups


class Hotspots:
    """Hotspots object represents a set of annotated hotspots.

    Attributes:
        hotspots (pd.Series): Series (list) of Hotspot objects.
        annotation (pd.DataFrame): Annotation table with description and statistics of hotspots.
        communities (dict): Dictionary representing annotated communities of hotspots
            (key - community_id, value - list of hotspot ids)
        island_rep_proteins_fasta (str): Path to a fasta file containing island representative proteins.
        prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

    """

    def __init__(self, hotspots: pd.Series, annotation: pd.DataFrame, parameters: ilund4u.manager.Parameters):
        """Hotspots class constructor.

        Arguments:
            hotspots (pd.Series): Series (list) of Hotspot objects.
            annotation (pd.DataFrame): Annotation table with description and statistics of hotspots.
            parameters (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        """
        self.hotspots = hotspots
        self.annotation = annotation
        self.communities = dict()
        self.prms = parameters
        self.__id_to_ind = {iid: idx for idx, iid in enumerate(self.annotation.index)}
        self.island_rep_proteins_fasta = os.path.join(parameters.args["output_dir"], "island_rep_proteins.fa")

    def save_as_db(self, db_folder: str) -> None:
        """Save Hotspots to the iLnd4u database.

        Arguments:
            db_folder (str): Database folder path.

        Returns:
            None

        """
        try:
            attributes_to_ignore = ["hotspots", "annotation", "communities_annot", "prms"]
            attributes = {k: v for k, v in self.__dict__.items() if k not in attributes_to_ignore}
            attributes["island_rep_proteins_fasta"] = os.path.basename(attributes["island_rep_proteins_fasta"])
            with open(os.path.join(db_folder, "hotspots.attributes.json"), 'w') as json_file:
                json.dump(attributes, json_file)
            self.annotation.to_csv(os.path.join(db_folder, "hotspots.annotations.tsv"), sep="\t",
                                   index_label="hotspot_id")
            island_annotation_table = pd.DataFrame()
            hotspot_db_ind = []
            for hotspot in self.hotspots.to_list():
                hotspot_db_ind.append(hotspot.get_hotspot_db_row())
                h_island_annot = hotspot.island_annotation.copy()
                h_island_annot["hotspot_id"] = hotspot.hotspot_id
                island_annotation_table = pd.concat([island_annotation_table, h_island_annot])

            with open(os.path.join(db_folder, "hotspot.ind.attributes.json"), "w") as json_file:
                json.dump(hotspot_db_ind, json_file)

            island_annotation_table.to_csv(os.path.join(db_folder, "hotspot.ind.island.annotations.tsv"), sep="\t",
                                           index_label="island")
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to write hotspots to the database.") from error

    @classmethod
    def db_init(cls, db_path: str, proteomes: Proteomes, parameters: ilund4u.manager.Parameters):
        """Class method to load a Proteomes object from a database.

        Arguments:
            db_path (str): path to the database.
            proteomes (Proteomes): Proteomes object.
            parameters (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        Returns:
            cls: Hotspots object.

        """
        try:
            if parameters.args["verbose"]:
                print(f"○ Loading hotspot objects...", file=sys.stdout)
            island_annotation = pd.read_table(os.path.join(db_path, "hotspot.ind.island.annotations.tsv"),
                                              sep="\t").set_index("island")
            with open(os.path.join(db_path, "hotspot.ind.attributes.json"), "r") as json_file:
                hotspot_ind_attributes = json.load(json_file)
            if parameters.args["verbose"]:
                bar = progress.bar.FillingCirclesBar(" ", max=len(hotspot_ind_attributes), suffix='%(index)d/%(max)d')
            hotspot_list = []
            for hotspot_dict in hotspot_ind_attributes:
                if parameters.args["verbose"]:
                    bar.next()
                hotspot_dict["island_annotation"] = island_annotation[
                    island_annotation["hotspot_id"] == hotspot_dict["hotspot_id"]].copy()
                hotspot_proteomes = proteomes.proteomes.loc[
                    proteomes.communities[hotspot_dict["proteome_community"]]].to_list()
                islands_list = [island for proteome in hotspot_proteomes for island in proteome.islands.to_list()]
                islands_series = pd.Series(islands_list, index=[island.island_id for island in islands_list])
                hotspot_dict["islands"] = islands_series.loc[hotspot_dict["island_annotation"].index].to_list()
                hotspot_list.append(Hotspot(**hotspot_dict))
            if parameters.args["verbose"]:
                bar.finish()
            hotspots = pd.Series(hotspot_list, index=[hotspot.hotspot_id for hotspot in hotspot_list])
            annotation = pd.read_table(os.path.join(db_path, "hotspots.annotations.tsv"),
                                       sep="\t", dtype={"community": "Int32"}).set_index("hotspot_id")
            cls_obj = cls(hotspots, annotation, parameters)
            with open(os.path.join(db_path, "hotspots.attributes.json"), "r") as json_file:
                attributes = json.load(json_file)
            cls_obj.communities = {int(k): v for k, v in attributes["communities"].items()}
            cls_obj.island_rep_proteins_fasta = os.path.join(db_path, attributes["island_rep_proteins_fasta"])
            return cls_obj
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to read hotspots from the database.") from error

    def pyhmmer_annotation(self, proteomes: Proteomes) -> None:
        """Run pyhhmmer hmmscan against a set of databases for additional annotation of hotspot proteins.

        Arguments:
            proteomes (Proteomes): Proteomes object.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Preparing data for additional island protein annotation with pyhmmer hmmscan...",
                      file=sys.stdout)
            hotspots_repr_proteins = set()
            for hotspot in self.hotspots.to_list():
                for island in hotspot.islands:
                    proteome = proteomes.proteomes.at[island.proteome]
                    isl_groups = island.get_locus_groups(proteome.cdss)
                    hotspots_repr_proteins.update(isl_groups)
            initial_fasta_file = Bio.SeqIO.index(proteomes.proteins_fasta_file, "fasta")
            with open(self.island_rep_proteins_fasta, "wb") as out_handle:
                for acc in hotspots_repr_proteins:
                    try:
                        out_handle.write(initial_fasta_file.get_raw(acc))
                    except:
                        pass

            alignment_table = ilund4u.methods.run_pyhmmer(self.island_rep_proteins_fasta, len(hotspots_repr_proteins),
                                                          self.prms)
            if not alignment_table.empty:
                found_hits_for = alignment_table.index.to_list()
                for proteome in proteomes.proteomes.to_list():
                    proteome_cdss = proteome.cdss.to_list()
                    proteome_cdss_with_hits = [cds.cds_id for cds in proteome_cdss if cds.group in found_hits_for]
                    if proteome_cdss_with_hits:
                        cdss_with_hits = proteome.cdss.loc[proteome_cdss_with_hits].to_list()
                        for cds in cdss_with_hits:
                            alignment_table_row = alignment_table.loc[cds.group]
                            cds.hmmscan_results = dict(db=alignment_table_row["target_db"],
                                                       target=alignment_table_row["t_description"],
                                                       evalue=alignment_table_row["hit_evalue"])
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to run pyhmmer hmmscan annotation.") from error

    def build_hotspot_network(self):
        """Build hotspot network and merge similar hotspots from different proteome communities into hotspot communities.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Hotspot network construction...", file=sys.stdout)
            hotspot_signature_sizes = pd.Series()
            hotspot_proteome_community = dict()
            signature_cluster_to_hotspot = collections.defaultdict(collections.deque)
            for hid, hotspot in enumerate(self.hotspots):
                hotspot_proteome_community[hid] = hotspot.proteome_community
                hotspot_signature_sizes.at[hid] = len(hotspot.conserved_signature)
                for cs_cluster in hotspot.conserved_signature:
                    signature_cluster_to_hotspot[cs_cluster].append(hid)
            edges, weights = [], []
            bar = progress.bar.FillingCirclesBar(" ", max=len(self.hotspots), suffix="%(index)d/%(max)d")
            for i, hotspot_i in enumerate(self.hotspots):
                bar.next()
                pc_i = hotspot_i.proteome_community
                signature_i = hotspot_i.conserved_signature
                size_i = hotspot_signature_sizes.iat[i]
                counts_i = collections.defaultdict(int)
                for sc in signature_i:
                    js = signature_cluster_to_hotspot[sc]
                    for j in js.copy():
                        if i < j:
                            if pc_i != hotspot_proteome_community[j]:  # to think about it
                                counts_i[j] += 1
                        else:
                            js.popleft()
                weights_i = pd.Series(counts_i)
                connected_n_sizes = hotspot_signature_sizes.iloc[weights_i.index]
                norm_factor_i = pd.Series(0.5 * (size_i + connected_n_sizes) / (size_i * connected_n_sizes), \
                                          index=weights_i.index)
                weights_i = weights_i.mul(norm_factor_i)
                weights_i = weights_i[weights_i >= self.prms.args["hotspot_similarity_cutoff"]]
                for j, w in weights_i.items():
                    edges.append([i, j])
                    weights.append(round(w, 4))
            bar.finish()
            print("○ Hotspot network partitioning using the Leiden algorithm...")
            graph = igraph.Graph(len(hotspot_signature_sizes.index), edges, directed=False)
            graph.vs["index"] = hotspot_signature_sizes.index.to_list()
            graph.vs["hotspot_id"] = [h.hotspot_id for h in self.hotspots]
            graph.es["weight"] = weights
            graph.save(os.path.join(self.prms.args["output_dir"], f"hotspot_network.gml"))
            partition_leiden = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition,
                                                        resolution_parameter=self.prms.args[
                                                            "leiden_resolution_parameter_h"],
                                                        weights="weight", n_iterations=-1)
            graph.vs["communities_Leiden"] = partition_leiden.membership
            hotspot_communities_annot_rows = []
            communities_sizes = []
            self.annotation["hotspot_community"] = pd.Series(dtype='Int64')
            for community_index, community in enumerate(partition_leiden):
                community_size = len(community)
                subgraph = graph.subgraph(community)
                hotspots = subgraph.vs["hotspot_id"]
                self.communities[community_index] = hotspots
                self.annotation.loc[hotspots, "hotspot_community"] = community_index
                if community_size > 1:
                    communities_sizes.append(community_size)
                    subgraph_edges = subgraph.get_edgelist()
                    num_of_edges = len(subgraph_edges)
                    num_of_edges_fr = num_of_edges / (community_size * (community_size - 1) * 0.5)
                    weights = subgraph.es["weight"]
                    avg_weight = round(np.mean(weights), 3)
                    max_identity = max(weights)
                else:
                    num_of_edges, num_of_edges_fr, avg_weight, max_identity = "", "", "", ""
                hotspot_communities_annot_rows.append([community_index, community_size, avg_weight, max_identity,
                                                       num_of_edges_fr, ";".join(hotspots)])
            communities_annot = pd.DataFrame(hotspot_communities_annot_rows, columns=["id", "size", "avg_weight",
                                                                                      "max_weight", "fr_edges",
                                                                                      "hotspots"])
            communities_annot.to_csv(os.path.join(os.path.join(self.prms.args["output_dir"],
                                                               "hotspot_communities.tsv")),
                                     sep="\t", index=False)
            if self.prms.args["verbose"]:
                print(f"  ⦿  {sum(communities_sizes)} hotspots were merged to {len(communities_sizes)} not singleton "
                      f"communities")
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to build hotspot network.") from error

    def calculate_hotspot_statistics_and_get_annotation(self, proteomes: Proteomes) -> pd.DataFrame:
        """Calculate hotspot statistics based using hmmscan results and save annotation tables.

        Arguments:
            proteomes (Proteomes): Proteomes object.

        Returns:
            pd.DataFrame: hotspot community annotation table.

        """
        try:
            hotspot_community_annot_rows = []
            r_types = ["cargo", "flanking"]
            # Create new columns
            for r_type in r_types:
                self.annotation[f"N_{r_type}_groups"] = pd.Series(dtype='Int64')
                for db_name in self.prms.args["database_names"]:
                    self.annotation[f"N_{db_name}_{r_type}_groups"] = pd.Series(dtype='Int64')
            # Get stat
            for h_com, hotspot_ids in self.communities.items():
                h_com_stat = collections.defaultdict(lambda: collections.defaultdict(dict))
                hotspot_com_groups = dict(cargo=set(), flanking=set())
                hotspots = self.hotspots.loc[hotspot_ids].to_list()
                n_islands = 0
                for hotspot in hotspots:
                    n_islands += hotspot.size
                    hotspot_groups = hotspot.get_hotspot_groups(proteomes)
                    db_stat = hotspot.calculate_database_hits_stats(proteomes, self.prms)
                    for r_type in r_types:
                        hotspot_com_groups[r_type].update(hotspot_groups[r_type])
                        self.annotation.at[hotspot.hotspot_id, f"N_{r_type}_groups"] = len(hotspot_groups[r_type])
                    for db_name in self.prms.args["database_names"]:
                        for r_type in r_types:
                            h_com_stat[db_name][r_type].update(db_stat[db_name][r_type])
                            self.annotation.at[hotspot.hotspot_id, f"N_{db_name}_{r_type}_groups"] = \
                                len(set(db_stat[db_name][r_type].values()))
                hc_annot_row = dict(com_id=h_com, community_size=len(hotspot_ids), N_islands=n_islands,
                                    hotspots=",".join(hotspot_ids), pdf_filename=f"{'_'.join(hotspot_ids)}.pdf")
                for r_type in r_types:
                    hc_annot_row[f"N_{r_type}_groups"] = len(hotspot_com_groups[r_type])
                for db_name in self.prms.args["database_names"]:
                    for r_type in r_types:
                        hc_annot_row[f"N_{db_name}_{r_type}_groups"] = len(set(h_com_stat[db_name][r_type].values()))
                hotspot_community_annot_rows.append(hc_annot_row)

            hotspot_community_annot = pd.DataFrame(hotspot_community_annot_rows)
            for db_name in self.prms.args["database_names"]:
                self.annotation[f"{db_name}_cargo_normalised"] = \
                    self.annotation.apply(lambda row: round(row[f"N_{db_name}_cargo_groups"] / row[f"N_cargo_groups"],
                                                            4), axis=1)
                hotspot_community_annot[f"{db_name}_cargo_normalised"] = \
                    hotspot_community_annot.apply(
                        lambda row: round(row[f"N_{db_name}_cargo_groups"] / row[f"N_cargo_groups"], 4), axis=1)
            self.annotation.to_csv(os.path.join(self.prms.args["output_dir"], "hotspot_annotation.tsv"), sep="\t",
                                   index_label="hotspot_id")
            hotspot_community_annot.to_csv(os.path.join(self.prms.args["output_dir"],
                                                        "hotspot_community_annotation.tsv"), sep="\t", index=False)

            # Save island annotation table
            island_annotation_table = pd.DataFrame()
            for hotspot in self.hotspots.to_list():
                h_island_annot = hotspot.island_annotation.copy()
                h_island_annot = h_island_annot.drop(columns=["degree", "strength", "island_index"])
                h_island_annot["hotspot_id"] = hotspot.hotspot_id
                columns = list(h_island_annot.columns)
                columns.insert(0, columns.pop(columns.index("hotspot_id")))
                h_island_annot = h_island_annot[columns]
                island_annotation_table = pd.concat([island_annotation_table, h_island_annot])

            island_annotation_table.to_csv(os.path.join(self.prms.args["output_dir"], "island_annotation.tsv"),
                                           sep="\t", index_label="island")

            return hotspot_community_annot
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to calculate hotspot and hotspot community statistics based "
                                               "on hmmscan results") from error

    def get_each_protein_group_statistics(self, proteomes: Proteomes) -> pd.DataFrame:
        """Calculate statistics for each protein group that includes its "jumping" properties and types of hotspots
            where it's encoded.

        Arguments:
            proteomes (Proteomes): Proteomes object.

        Returns:
            pd.DataFrame: Protein group statistics table

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Protein group statistics calculation...", file=sys.stdout)
                bar = progress.bar.FillingCirclesBar(" ", max=len(self.hotspots.index), suffix='%(index)d/%(max)d')
            protein_group_statistics_dict = collections.defaultdict(
                lambda: {"Hotspot_communities": set(), "Hotspots": set(), "Islands": set(),
                         "Proteome_communities": set(),
                         "Counts": 0, "db": "None", "db_hit": "None", "Name": "", "RepLength": 0})
            for h_com, hotspots_ids in self.communities.items():
                hotspots = self.hotspots.loc[hotspots_ids].to_list()
                for hotspot in hotspots:
                    if self.prms.args["verbose"]:
                        bar.next()
                    # Hotspot stat
                    for island in hotspot.islands:
                        proteome = proteomes.proteomes.at[island.proteome]
                        island_proteins = island.get_island_proteins(proteome.cdss)
                        island_protein_groups = island.get_island_groups(proteome.cdss)
                        for isp, ispg in zip(island_proteins, island_protein_groups):
                            if protein_group_statistics_dict[ispg]["Counts"] == 0:
                                cds_obj = proteome.cdss.at[isp]
                                protein_group_statistics_dict[ispg]["Name"] = cds_obj.name
                                protein_group_statistics_dict[ispg]["RepLength"] = cds_obj.length
                                if cds_obj.hmmscan_results:
                                    protein_group_statistics_dict[ispg]["db_hit"] = cds_obj.hmmscan_results["target"]
                                    protein_group_statistics_dict[ispg]["db"] = cds_obj.hmmscan_results["db"]
                            protein_group_statistics_dict[ispg]["Counts"] += 1
                            protein_group_statistics_dict[ispg]["Hotspot_communities"].add(h_com)
                            protein_group_statistics_dict[ispg]["Proteome_communities"].add(hotspot.proteome_community)
                            protein_group_statistics_dict[ispg]["Hotspots"].add(hotspot.hotspot_id)
                            protein_group_statistics_dict[ispg]["Islands"].add(island.island_id)
            if self.prms.args["verbose"]:
                bar.finish()
            statistic_rows = []
            for pg, pg_dict in protein_group_statistics_dict.items():
                row_dict = dict(representative_protein=pg, db=pg_dict["db"], pb_hit=pg_dict["db_hit"],
                                name=pg_dict["Name"], counts=pg_dict["Counts"], length=pg_dict["RepLength"])
                for k, v in pg_dict.items():
                    if isinstance(v, set):
                        row_dict[f"N_{k}"] = len(v)
                for db_name in self.prms.args["database_names"]:
                    dbsn_norm_values = []
                    for hotspot in pg_dict["Hotspots"]:
                        dbsn_norm_value = self.annotation.at[hotspot, f"{db_name}_cargo_normalised"]
                        if pg_dict["db"] == db_name:
                            dbsn_norm_value -= 1 / self.annotation.at[hotspot, f"N_cargo_groups"]
                        dbsn_norm_values.append(dbsn_norm_value)
                    row_dict[f"{db_name}_avg_cargo_fraction"] = round(np.mean(dbsn_norm_values), 4)
                    row_dict[f"{db_name}_max_cargo_fraction"] = round(max(dbsn_norm_values), 4)
                row_dict["hotspots"] = ",".join(pg_dict["Hotspots"])

                statistic_rows.append(row_dict)
            statistic_table = pd.DataFrame(statistic_rows)
            statistic_table.sort_values(by=["N_Hotspot_communities", "N_Hotspots"], ascending=[False, False],
                                        inplace=True)
            statistic_table.to_csv(os.path.join(self.prms.args["output_dir"], "protein_group_statistics.tsv"),
                                   sep="\t", index=False)
            return statistic_table
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to calculate protein group statistics.") from error


class Database:
    """Database object represents iLund4u database with proteomes and hotspots objects.

    Attributes:
        proteomes (Proteomes): Database proteomes object.
        hotspots (Hotspots): Database hotspots object.
        db_paths (dict): Dictionary of database paths.
        prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

    """

    def __init__(self, proteomes: Proteomes, hotspots: Hotspots, db_paths: dict,
                 parameters: ilund4u.manager.Parameters):
        """Database class constructor.

        Args:
            proteomes (Proteomes): Database proteomes object.
            hotspots (Hotspots): Database hotspots object.
            db_paths (dict): Dictionary of database paths.
            prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        """
        self.proteomes = proteomes
        self.hotspots = hotspots
        self.db_paths = db_paths
        self.prms = parameters

    def mmseqs_search_versus_protein_database(self, query_fasta: str) -> pd.DataFrame:
        """Run mmseqs search versus protein database.

        Arguments:
            query_fasta (str): path to a query fasta file with protein sequence(s).

        Returns:
            pd.DataFrame: mmseqs search results table.

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Running mmseqs for protein search versus the database of proteins...",
                      file=sys.stdout)
            if not os.path.exists(self.prms.args["output_dir"]):
                os.mkdir(self.prms.args["output_dir"])
            mmseqs_output_folder = os.path.join(self.prms.args["output_dir"], "mmseqs")
            if os.path.exists(mmseqs_output_folder):
                shutil.rmtree(mmseqs_output_folder)
            os.mkdir(mmseqs_output_folder)
            mmseqs_output_folder_db = os.path.join(mmseqs_output_folder, "DBs")
            os.mkdir(mmseqs_output_folder_db)
            mmseqs_stdout = open(os.path.join(mmseqs_output_folder, "mmseqs_stdout.txt"), "w")
            mmseqs_stderr = open(os.path.join(mmseqs_output_folder, "mmseqs_stderr.txt"), "w")
            query_length = len(list(Bio.SeqIO.parse(query_fasta, "fasta")))
            subprocess.run([self.prms.args["mmseqs_binary"], "createdb", query_fasta,
                            os.path.join(mmseqs_output_folder_db, "query_seq_db")], stdout=mmseqs_stdout,
                           stderr=mmseqs_stderr)
            subprocess.run([self.prms.args["mmseqs_binary"], "search",
                            os.path.join(mmseqs_output_folder_db, "query_seq_db"),
                            self.db_paths["proteins_db"],
                            os.path.join(mmseqs_output_folder_db, "search_res_db"),
                            os.path.join(mmseqs_output_folder, "tmp"), "-e",
                            str(self.prms.args["mmseqs_search_evalue"]),
                            "-s",  str(self.prms.args["mmseqs_search_s"])], stdout=mmseqs_stdout, stderr=mmseqs_stderr)
            subprocess.run([self.prms.args["mmseqs_binary"], "convertalis",
                            os.path.join(mmseqs_output_folder_db, "query_seq_db"),
                            self.db_paths["proteins_db"],
                            os.path.join(mmseqs_output_folder_db, "search_res_db"),
                            os.path.join(mmseqs_output_folder, "mmseqs_search_results.tsv"), "--format-output",
                            "query,target,qlen,tlen,alnlen,fident,qstart,qend,tstart,tend,evalue",
                            "--format-mode", "4"], stdout=mmseqs_stdout, stderr=mmseqs_stderr)
            mmseqs_search_results = pd.read_table(os.path.join(mmseqs_output_folder, "mmseqs_search_results.tsv"),
                                                  sep="\t")
            mmseqs_search_results["qcov"] = mmseqs_search_results.apply(lambda row: row["alnlen"] / row["qlen"], axis=1)
            mmseqs_search_results["tcov"] = mmseqs_search_results.apply(lambda row: row["alnlen"] / row["tlen"], axis=1)
            mmseqs_search_results = mmseqs_search_results[
                (mmseqs_search_results["qcov"] >= self.prms.args["mmseqs_search_qcov"]) &
                (mmseqs_search_results["tcov"] >= self.prms.args["mmseqs_search_tcov"]) &
                (mmseqs_search_results["fident"] >= self.prms.args["mmseqs_search_fident"])]
            queries_with_res = len(set(mmseqs_search_results["query"].to_list()))
            target_to_group = dict()
            for proteome in self.proteomes.proteomes.to_list():
                for cds in proteome.cdss.to_list():
                    target_to_group[cds.cds_id] = cds.group
            mmseqs_search_results["group"] = mmseqs_search_results["target"].apply(lambda t: target_to_group[t])
            mmseqs_search_results.to_csv(os.path.join(self.prms.args["output_dir"], "mmseqs_homology_search_full.tsv"),
                                         sep="\t", index=False)
            if self.prms.args["verbose"]:
                if queries_with_res > 0:
                    print(f"  ⦿ A homologous group was found for {queries_with_res}/{query_length} query protein"
                          f"{'s' if query_length > 1 else ''}", file=sys.stdout)
                else:
                    print(f"  ⦿ No homologous group was found for {query_length} query protein"
                          f"{'s' if query_length > 1 else ''}", file=sys.stdout)
            return mmseqs_search_results
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to run mmseqs search versus protein database.") from error

    def protein_search_mode(self, query_fasta: str, query_label: typing.Union[None, str] = None) -> None:
        """Run protein search mode which finds homologues of your query proteins in the database and returns
            comprehensive output including visualisation and hotspot annotation.

        Arguments:
            query_fasta (str): Fasta with query protein sequence.
            query_label (str): Label to be shown on lovis4u visualisation.

        Returns:
            None
        """
        try:
            # Load fasta
            query_records = list(Bio.SeqIO.parse(query_fasta, "fasta"))
            if len(query_records) > 1:
                raise ilund4u.manager.ilund4uError("Only single query protein is allowed for protein mode")
            query_record = query_records[0]
            # Run mmseqs for homology search
            mmseqs_results = self.mmseqs_search_versus_protein_database(query_fasta)
            if len(mmseqs_results.index) == 0:
                print("○ Termination since no homology to hotspot db proteins was found", file=sys.stdout)
                return None
            if self.prms.args["protein_search_target_mode"] == "proteins":
                homologous_protein_ids = mmseqs_results["target"].to_list()
            elif self.prms.args["protein_search_target_mode"] == "group":
                mmseqs_results.sort_values(by=["evalue", "qcov", "tcov", "fident"],
                                           ascending=[True, False, False, False], inplace=True)
                mmseqs_results = mmseqs_results.drop_duplicates(subset="query", keep="first").set_index("query")
                homologous_protein_ids = []
                homologous_group = mmseqs_results.at[query_record.id, "group"]
            # Searching for hotspots
            if self.prms.args["verbose"]:
                print(f"○ Searching for hotspots with your query protein homologues...", file=sys.stdout)
            found_hotspots = collections.defaultdict(list)
            island_annotations = []
            location_stat = dict(flanking=0, cargo=0)
            for hotspot in self.hotspots.hotspots.to_list():
                for island in hotspot.islands:
                    proteome = self.proteomes.proteomes.at[island.proteome]
                    if self.prms.args["protein_search_target_mode"] == "proteins":
                        locus_proteins = island.get_locus_proteins(proteome.cdss)
                        overlapping = list(set(homologous_protein_ids) & set(locus_proteins))
                    elif self.prms.args["protein_search_target_mode"] == "group":
                        locus_groups = island.get_locus_groups(proteome.cdss)
                        overlapping = homologous_group in locus_groups
                    if overlapping:
                        if self.prms.args["protein_search_target_mode"] == "group":
                            homologous_protein_ids_island = []
                            locus_proteins = island.get_locus_proteins(proteome.cdss)
                            for lp, lpg in zip(locus_proteins, locus_groups):
                                if lpg == homologous_group:
                                    homologous_protein_ids_island.append(lp)
                                    homologous_protein_ids.append(lp)
                            overlapping = homologous_protein_ids_island

                        isl_groups = island.get_island_groups(proteome.cdss)
                        isl_proteins = island.get_island_proteins(proteome.cdss)
                        for op in overlapping:
                            if op in isl_proteins:
                                location_stat["cargo"] += 1
                            else:
                                location_stat["flanking"] += 1
                        island_annotation = hotspot.island_annotation.loc[island.island_id].copy()
                        island_annotation.drop(labels=["island_index", "strength", "degree"], inplace=True)
                        island_annotation.at["indexes"] = ",".join(map(str, island.indexes))
                        island_annotation.at["size"] = island.size
                        island_annotation.at["island_proteins"] = ",".join(island.get_island_proteins(proteome.cdss))
                        island_annotation.at["island_protein_groups"] = ",".join(isl_groups)
                        island_annotation.at["query_homologues"] = ",".join(overlapping)
                        island_annotations.append(island_annotation)
                        found_hotspots[hotspot.hotspot_id].append(island)
            if sum(location_stat.values()) == 0:
                print("○ Termination since no homologous protein was found in hotspots (neither as flanking or cargo)",
                      file=sys.stdout)
                return None
            found_islands = [island.island_id for islands in found_hotspots.values() for island in islands]
            island_annotations = pd.DataFrame(island_annotations)
            island_annotations.to_csv(os.path.join(os.path.join(self.prms.args["output_dir"],
                                                                "found_island_annotation.tsv")),
                                      sep="\t", index_label="island_id")
            found_hotspots_annotation = self.hotspots.annotation.loc[found_hotspots.keys()]
            found_hotspots_annotation.to_csv(os.path.join(self.prms.args["output_dir"], "found_hotspot_annotation.tsv"),
                                             sep="\t", index_label="hotspot_id")
            found_hotspot_communities = list(set(found_hotspots_annotation["hotspot_community"].to_list()))
            if self.prms.args["verbose"]:
                print(f"  ⦿ Query protein homologues were found in {len(found_hotspot_communities)} hotspot "
                      f"communit{'y' if len(found_hotspot_communities) == 1 else 'ies'} "
                      f"({len(found_hotspots.keys())} hotspot{'s' if len(found_hotspots.keys()) > 1 else ''}) on "
                      f"{len(found_islands)} island{'s' if len(found_islands) > 1 else ''}\n"
                      f"    Found as cargo: {location_stat['cargo']}, as flanking gene: {location_stat['flanking']}",
                      file=sys.stdout)
            homologous_protein_fasta = os.path.join(self.prms.args["output_dir"], "homologous_proteins.fa")
            full_fasta_file = Bio.SeqIO.index(self.proteomes.proteins_fasta_file, "fasta")
            with open(homologous_protein_fasta, "w") as out_handle:
                Bio.SeqIO.write(query_record, out_handle, "fasta")
                for acc in homologous_protein_ids:
                    out_handle.write(full_fasta_file.get_raw(acc).decode())
            # MSA visualisation
            msa4u_p = msa4u.manager.Parameters()
            msa4u_p.arguments["label"] = "id"
            msa4u_p.arguments["verbose"] = False
            msa4u_p.arguments["output_filename"] = os.path.join(self.prms.args["output_dir"],
                                                                "msa4u_homologous_proteines.pdf")
            msa4u_p.arguments["output_filename_aln"] = os.path.join(self.prms.args["output_dir"],
                                                                    "homologous_proteins_aln.fa")
            fasta = msa4u.manager.Fasta(fasta=homologous_protein_fasta, parameters=msa4u_p)
            mafft_output = fasta.run_mafft()
            msa = msa4u.manager.MSA(mafft_output, msa4u_p)
            msa.plot()
            if self.prms.args["verbose"]:
                print(f"⦿ Homologous proteins were saved to {homologous_protein_fasta} and the MSA was visualised with "
                      f"MSA4u")
                print(f"○ Visualisation of the hotspot(s) with your query protein homologues using lovis4u...",
                      file=sys.stdout)
            # lovis4u visualisation
            vis_output_folders = [os.path.join(self.prms.args["output_dir"], "lovis4u_full"),
                                  os.path.join(self.prms.args["output_dir"], "lovis4u_with_query")]
            for vis_output_folder in vis_output_folders:
                if os.path.exists(vis_output_folder):
                    shutil.rmtree(vis_output_folder)
                os.mkdir(vis_output_folder)
            additional_annotation = dict()
            for hpid in homologous_protein_ids:
                additional_annotation[hpid] = dict(stroke_colour="#000000", fill_colour="#000000")
                if query_label:
                    additional_annotation[hpid]["name"] = query_label
            drawing_manager = ilund4u.drawing.DrawingManager(self.proteomes, self.hotspots, self.prms)
            for community in found_hotspot_communities:
                drawing_manager.plot_hotspots(self.hotspots.communities[community],
                                              output_folder=os.path.join(self.prms.args["output_dir"], "lovis4u_full"),
                                              additional_annotation=additional_annotation)
            drawing_manager.plot_hotspots(found_hotspots.keys(),
                                          output_folder=os.path.join(self.prms.args["output_dir"],
                                                                     "lovis4u_with_query"),
                                          island_ids=found_islands,
                                          additional_annotation=additional_annotation)
            if self.prms.args["verbose"]:
                print(f"⦿ Done!")
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to perform protein search versus the database.") from error

    def proteome_annotation_mode(self, query_gff: str) -> None:
        """Run proteome annotation mode which searches for similar proteomes in the database and annotate hotspots and
            variable proteins in the query proteome in case a community with similar proteomes was found in the database.

        Arguments:
            query_gff (str): GFF with query proteome.

        Returns:
            None
        """
        try:
            # Load query gff
            proteomes_helper_obj = ilund4u.data_processing.Proteomes(parameters=self.prms)
            proteomes_helper_obj.load_sequences_from_extended_gff(input_f=[query_gff])
            query_proteome = proteomes_helper_obj.proteomes.iat[0]
            # Get and parse mmseqs search results
            mmseqs_results = self.mmseqs_search_versus_protein_database(proteomes_helper_obj.proteins_fasta_file)
            if self.prms.args["verbose"]:
                print(f"○ Searching for similar proteomes in the database network", file=sys.stdout)
            mmseqs_results.sort_values(by=["evalue", "qcov", "tcov", "fident"], ascending=[True, False, False, False],
                                       inplace=True)
            mmseqs_results = mmseqs_results.drop_duplicates(subset="query", keep="first").set_index("query")
            proteins_wo_hits = []
            for cds in query_proteome.cdss.to_list():
                if cds.cds_id in mmseqs_results.index:
                    cds.group = mmseqs_results.at[cds.cds_id, "group"]
                else:
                    cds.group = f"{cds.cds_id}"
                    proteins_wo_hits.append(cds.group)
            # Running pyhmmer annotation
            if self.prms.args["verbose"]:
                print(f"○ Preparing data for protein annotation with pyhmmer hmmscan...", file=sys.stdout)
            alignment_table = ilund4u.methods.run_pyhmmer(proteomes_helper_obj.proteins_fasta_file,
                                                          len(query_proteome.cdss.index), self.prms)
            if not alignment_table.empty:
                found_hits_for = alignment_table.index.to_list()
                proteome_cdss = query_proteome.cdss.to_list()
                proteome_cdss_with_hits = [cds.cds_id for cds in proteome_cdss if cds.cds_id in found_hits_for]
                if proteome_cdss_with_hits:
                    cdss_with_hits = query_proteome.cdss.loc[proteome_cdss_with_hits].to_list()
                    for cds in cdss_with_hits:
                        alignment_table_row = alignment_table.loc[cds.cds_id]
                        cds.hmmscan_results = dict(db=alignment_table_row["target_db"],
                                                   target=alignment_table_row["t_description"],
                                                   evalue=alignment_table_row["hit_evalue"])
            # Connect to the database proteome network
            proteome_names = pd.Series({idx: sid for idx, sid in enumerate(self.proteomes.annotation.index)})
            proteome_sizes = self.proteomes.annotation[["proteome_size_unique", "index"]]
            proteome_sizes = proteome_sizes.set_index("index")["proteome_size_unique"]
            cluster_to_sequences = collections.defaultdict(list)
            for p_index, proteome in enumerate(self.proteomes.proteomes.to_list()):
                cds_groups = set(proteome.cdss.apply(lambda cds: cds.group).to_list())
                for cds_g in cds_groups:
                    cluster_to_sequences[cds_g].append(proteome.proteome_id)
            cluster_to_proteome_index = dict()
            for cluster, sequences in cluster_to_sequences.items():
                indexes = sorted([self.proteomes.seq_to_ind[seq_id] for seq_id in sequences])
                cluster_to_proteome_index[cluster] = indexes
            query_clusters = set(query_proteome.cdss.apply(lambda cds: cds.group).to_list())
            query_size = len(query_clusters)
            counts = collections.defaultdict(int)
            for cl in query_clusters:
                if cl not in proteins_wo_hits:
                    js = cluster_to_proteome_index[cl]
                    for j in js:
                        counts[j] += 1
            weights = pd.Series(counts)
            proteome_sizes_connected = proteome_sizes.iloc[weights.index]
            norm_factor = pd.Series(
                0.5 * (query_size + proteome_sizes_connected) / (query_size * proteome_sizes_connected), \
                index=weights.index)
            weights = weights.mul(norm_factor)
            weights = weights[weights >= self.prms.args["proteome_similarity_cutoff"]]
            query_network_df = pd.DataFrame(dict(weight=weights, seq_id=proteome_names.iloc[weights.index]))
            query_network_df.sort_values(by="weight", inplace=True, ascending=False)
            query_network_df.to_csv(os.path.join(self.prms.args["output_dir"], "query_proteome_network.tsv"),
                                    sep="\t", index_label="t_index")
            query_network_df = query_network_df.set_index("seq_id")
            max_weight = round(query_network_df["weight"].max(), 2)
            if len(weights.index) == 0:
                print("○ Termination since no similar proteome was found in the database", file=sys.stdout)
                sys.exit()
            if self.prms.args["verbose"]:
                print(f"⦿ {len(weights.index)} similar proteomes were found in the database network with "
                      f"max proteome similarity = {max_weight}", file=sys.stdout)
            similar_proteoms_ids = query_network_df.index.to_list()
            # Assign the closest community
            similar_communities_rows = list()
            for pcom_id, pcom_pr_ids in self.proteomes.communities.items():
                com_size = len(pcom_pr_ids)
                overlapping = list(set(pcom_pr_ids) & set(similar_proteoms_ids))
                if overlapping:
                    weights_subset = query_network_df.loc[overlapping]["weight"].to_list()
                    similar_communities_rows.append(dict(com_id=pcom_id, com_size=com_size,
                                                         connection_fr=len(overlapping) / com_size,
                                                         avg_weight=np.mean(weights_subset)))
            similar_communities = pd.DataFrame(similar_communities_rows)
            similar_communities.sort_values(by=["avg_weight", "connection_fr", "com_size"], inplace=True,
                                            ascending=[False, False, False])
            similar_communities.to_csv(os.path.join(self.prms.args["output_dir"], "similar_proteome_communities.tsv"),
                                       sep="\t", index=False)
            selected_community_dict = similar_communities.iloc[0].to_dict()
            selected_community_dict["com_id"] = int(selected_community_dict["com_id"])
            query_proteome_community = selected_community_dict["com_id"]
            if self.prms.args["verbose"]:
                print(f"⦿ The query proteome was assigned to a community (id: {int(selected_community_dict['com_id'])})"
                      f" with connection to {round(selected_community_dict['connection_fr'], 2)} of its members with "
                      f"avg weight = {round(selected_community_dict['avg_weight'], 2)}", file=sys.stdout)
            # Define cds classes
            com_protein_classes = dict()
            for com_proteome in self.proteomes.communities[selected_community_dict["com_id"]]:
                for cds in self.proteomes.proteomes.at[com_proteome].cdss:
                    com_protein_classes[cds.group] = cds.g_class
            defined_classes_rows = []
            for cds in query_proteome.cdss:
                if cds.group in com_protein_classes.keys():
                    cds.g_class = com_protein_classes[cds.group]
                else:
                    cds.g_class = "variable"
                defined_classes_rows.append(dict(cds_id=cds.cds_id, cds_class=cds.g_class))
            defined_classes = pd.DataFrame(defined_classes_rows)
            defined_classes.to_csv(os.path.join(self.prms.args["output_dir"], "query_protein_clusters.tsv"),
                                   sep="\t", index=False)
            defined_classes_c = collections.Counter(defined_classes["cds_class"].to_list())
            if self.prms.args["verbose"]:
                print(f"⦿ Protein class distribution in query proteome: "
                      f"{', '.join(f'{v} {k}' for k, v in defined_classes_c.items())}"
                      , file=sys.stdout)
            # Annotate variable islands
            query_proteome.annotate_variable_islands(self.prms)
            if self.prms.args["verbose"]:
                print(f"⦿ {len(query_proteome.islands.index)} variable island"
                      f"{'s were' if len(query_proteome.islands.index) > 1 else ' was'} annotated in the "
                      f"query proteome", file=sys.stdout)
            # Connect query proteome islands to the island network
            community_hotspots = [h for h in self.hotspots.hotspots.to_list() if
                                  h.proteome_community == query_proteome_community]
            if community_hotspots:
                hotspots_islands = [island for hotspot in community_hotspots for island in hotspot.islands]
                island_id_to_index = {isl.island_id: ind for ind, isl in enumerate(hotspots_islands)}
                com_island_n_sizes = pd.Series()
                com_neighbours = pd.Series()
                cluster_to_island = collections.defaultdict(list)
                for island in hotspots_islands:
                    island_proteome = self.proteomes.proteomes.at[island.proteome]
                    island_id = island.island_id
                    island_index = island_id_to_index[island_id]
                    conserved_island_neighbours_groups = set(island.get_cons_neighbours_groups(island_proteome.cdss))
                    com_neighbours.at[island_index] = list(conserved_island_neighbours_groups)
                    com_island_n_sizes.at[island_index] = len(conserved_island_neighbours_groups)
                    for cing in conserved_island_neighbours_groups:
                        cluster_to_island[cing].append(island_index)
                hotspot_hits_statistic_rows = []
                for q_island in query_proteome.islands.to_list():
                    q_isl_neighbours = list(set(q_island.get_cons_neighbours_groups(query_proteome.cdss)))
                    q_isl_size = len(q_isl_neighbours)
                    q_isl_counts = collections.defaultdict(int)
                    for qin in q_isl_neighbours:
                        js = cluster_to_island[qin]
                        for j in js:
                            q_isl_counts[j] += 1
                    q_isl_weights = pd.Series(q_isl_counts)
                    q_isl_connected_n_sizes = com_island_n_sizes.iloc[q_isl_weights.index]
                    q_isl_norm_factors = pd.Series(
                        0.5 * (q_isl_size + q_isl_connected_n_sizes) / (q_isl_size * q_isl_connected_n_sizes),
                        index=q_isl_weights.index)
                    q_isl_weights = q_isl_weights.mul(q_isl_norm_factors)
                    q_isl_weights = q_isl_weights[
                        q_isl_weights >= self.prms.args["island_neighbours_similarity_cutoff"]]
                    q_isl_similar_islands = q_isl_weights.index.to_list()
                    similar_hotspots_rows = list()
                    for c_hotspot in community_hotspots:
                        hotspot_islands = c_hotspot.islands
                        hotspot_islsnds_idx = [island_id_to_index[isl.island_id] for isl in hotspot_islands]
                        overlapping = list(set(q_isl_similar_islands) & set(hotspot_islsnds_idx))
                        if overlapping:
                            weights_subset = q_isl_weights.loc[overlapping].to_list()
                            connection_fr = len(overlapping) / c_hotspot.size
                            if connection_fr >= self.prms.args["island_neighbours_similarity_cutoff"]:
                                similar_hotspots_rows.append(dict(hotspot_id=c_hotspot.hotspot_id,
                                                                  hotspot_size=c_hotspot.size,
                                                                  connection_fr=connection_fr,
                                                                  avg_weight=np.mean(weights_subset)))
                    if similar_hotspots_rows:
                        similar_hotspots_rows = pd.DataFrame(similar_hotspots_rows)
                        similar_hotspots_rows.sort_values(by=["avg_weight", "connection_fr", "hotspot_size"],
                                                          inplace=True,
                                                          ascending=[False, False, False])
                        best_hit_for_q_isl = similar_hotspots_rows.iloc[0].to_dict()
                        hotspot_hits_statistic_rows.append(dict(query_island_id=q_island.island_id,
                                                                query_island_proteins=",".join(
                                                                    q_island.get_island_proteins(query_proteome.cdss)),
                                                                closest_hotspot=best_hit_for_q_isl["hotspot_id"],
                                                                closest_hotspot_size=best_hit_for_q_isl["hotspot_size"],
                                                                closest_hotspot_avg_weight=best_hit_for_q_isl[
                                                                    "avg_weight"],
                                                                closest_hotspot_connection_fr=best_hit_for_q_isl[
                                                                    "connection_fr"]))
                hotspot_hits_statistic = pd.DataFrame(hotspot_hits_statistic_rows)
                hotspot_hits_statistic.sort_values(by=["closest_hotspot_avg_weight", "closest_hotspot_connection_fr"],
                                                   inplace=True, ascending=[False, False])
                hotspot_hits_statistic.to_csv(
                    os.path.join(self.prms.args["output_dir"], "island_to_hotspot_mapping.tsv"),
                    sep="\t", index=False)
                hotspot_hits_statistic = hotspot_hits_statistic.drop_duplicates(subset="closest_hotspot", keep="first")
                hotspot_hits_statistic.to_csv(os.path.join(self.prms.args["output_dir"],
                                                           "island_to_hotspot_mapping_deduplicated.tsv"),
                                              sep="\t", index=False)
                subset_hotspot_annotation = self.hotspots.annotation.loc[hotspot_hits_statistic["closest_hotspot"]]
                subset_hotspot_annotation.to_csv(os.path.join(self.prms.args["output_dir"],
                                                              "annotation_of_mapped_hotspots.tsv"), sep="\t",
                                                 index=False)
                if self.prms.args["verbose"]:
                    print(f"⦿ {len(hotspot_hits_statistic.index)} of annotated variable island"
                          f"{'s were' if len(hotspot_hits_statistic.index) > 1 else ' was'} mapped to "
                          f"the database hotspots", file=sys.stdout)
            else:
                if self.prms.args["verbose"]:
                    print(f"⦿ The proteome community does not have any annotated hotspots", file=sys.stdout)
            # Update proteome and hotspots objects
            self.proteomes.proteomes.at[query_proteome.proteome_id] = query_proteome
            self.proteomes.annotation.loc[query_proteome.proteome_id] = proteomes_helper_obj.annotation.loc[
                query_proteome.proteome_id]
            self.proteomes.communities[query_proteome_community].append(query_proteome.proteome_id)
            if community_hotspots:
                for index, row in hotspot_hits_statistic.iterrows():
                    hotspot_id = row["closest_hotspot"]
                    query_island_id = row["query_island_id"]
                    query_island = query_proteome.islands.at[query_island_id]
                    closest_hotspot = self.hotspots.hotspots.at[hotspot_id]
                    closest_hotspot.islands.append(query_island)
            # Visualisation
            if self.prms.args["verbose"]:
                print(f"○ Lovis4u visualisation of communities and proteomes...", file=sys.stdout)
            drawing_manager = ilund4u.drawing.DrawingManager(self.proteomes, self.hotspots, self.prms)

            drawing_manager.plot_proteome_community(community=query_proteome_community,
                                                    output_folder=self.prms.args["output_dir"],
                                                    mode="hotspot",
                                                    filename="lovis4u_proteome_community_hotspots.pdf")
            drawing_manager.plot_proteome_community(community=query_proteome_community,
                                                    output_folder=self.prms.args["output_dir"],
                                                    mode="regular",
                                                    filename="lovis4u_proteome_community.pdf")
            drawing_manager.plot_proteome_community(community=query_proteome_community,
                                                    output_folder=self.prms.args["output_dir"],
                                                    mode="regular",
                                                    proteome_ids=[query_proteome.proteome_id],
                                                    filename="lovis4u_query_proteome_variable.pdf")
            drawing_manager.plot_proteome_community(community=query_proteome_community,
                                                    output_folder=self.prms.args["output_dir"],
                                                    mode="hotspot",
                                                    proteome_ids=[query_proteome.proteome_id],
                                                    filename="lovis4u_query_proteome_hotspot.pdf")
            if community_hotspots:
                if self.prms.args["verbose"]:
                    print(f"○ Lovis4u visualisation found hotspots..", file=sys.stdout)
                if len(hotspot_hits_statistic.index) > 0:
                    vis_output_folders = [os.path.join(self.prms.args["output_dir"], "lovis4u_hotspot_full"),
                                          os.path.join(self.prms.args["output_dir"], "lovis4u_hotspot_with_query")]
                    for vis_output_folder in vis_output_folders:
                        if os.path.exists(vis_output_folder):
                            shutil.rmtree(vis_output_folder)
                        os.mkdir(vis_output_folder)
                for index, row in hotspot_hits_statistic.iterrows():
                    hotspot_id = row["closest_hotspot"]
                    query_island = row["query_island_id"]
                    for hc, c_hotspots in self.hotspots.communities.items():
                        if hotspot_id in c_hotspots:
                            drawing_manager.plot_hotspots([hotspot_id],
                                                          output_folder=os.path.join(self.prms.args["output_dir"],
                                                                                     "lovis4u_hotspot_with_query"),
                                                          island_ids=[query_island])
                            drawing_manager.plot_hotspots(c_hotspots,
                                                          output_folder=os.path.join(self.prms.args["output_dir"],
                                                                                     "lovis4u_hotspot_full"),
                                                          keep_while_deduplication=[query_island])
            if self.prms.args["verbose"]:
                print(f"⦿ Done!")
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to run proteome annotation mode versus the database.") from error
