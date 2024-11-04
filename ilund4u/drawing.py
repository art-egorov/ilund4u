"""
This module provides visualisation classes and methods for the tool.
"""
import progress.bar
import pandas as pd
import tempfile
import typing
import sys
import os

import matplotlib.colors
import seaborn
import shutil

import lovis4u
import ilund4u.data_processing


class DrawingManager:
    """Manager for data visualisation using LoVis4u library.

    Attributes:
        proteomes (ilund4u.data_processing.Proteomes): Proteomes object.
        hotspots (ilund4u.data_processing.Hotspots): Hotspots object.
        prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.


    """

    def __init__(self, proteomes: ilund4u.data_processing.Proteomes, hotspots: ilund4u.data_processing.Hotspots,
                 parameters: ilund4u.manager.Parameters):
        """DrawingManager class constructor

        Arguments:
            proteomes (ilund4u.data_processing.Proteomes): Proteomes object.
            hotspots (ilund4u.data_processing.Hotspots): Hotspots object.
            parameters (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        """
        self.proteomes = proteomes
        self.hotspots = hotspots
        self.prms = parameters

    def plot_hotspot_communities(self, communities: typing.Union[None, list] = None,
                                 shortest_labels="auto") -> None:
        """Run visualisation of hotspot list for each hotspot community.

        Arguments:
            communities (None | list): list of communities to be plotted.
            shortest_labels (bool | auto): Whether to put 1-based cds index only instead of CDS id for hypothetical
                proteins.

        Returns:
            None

        """
        try:
            vis_output_folder = os.path.join(self.prms.args["output_dir"], "lovis4u_hotspots")
            if os.path.exists(vis_output_folder):
                shutil.rmtree(vis_output_folder)
            os.mkdir(vis_output_folder)
            if not communities:
                communities = self.hotspots.communities.values()
            if self.prms.args["verbose"]:
                print(f"○ Visualisation of hotspot communities using lovis4u...", file=sys.stdout)
            bar = progress.bar.FillingCirclesBar(" ", max=len(communities), suffix='%(index)d/%(max)d')
            for hc in communities:
                self.plot_hotspots(hc, vis_output_folder, shortest_labels=shortest_labels)
                bar.next()
            bar.finish()
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to plot hotspot communities.") from error

    def plot_hotspots(self, hotspot_ids: list, output_folder: str = "default",
                      island_ids: typing.Union[None, list] = None,
                      proteome_ids: typing.Union[None, list] = None,
                      additional_annotation: typing.Union[None, dict] = None, keep_while_deduplication: list = [],
                      shortest_labels: typing.Union[str, bool] = "auto", compact_mode: bool = False,
                      keep_temp_data=True):
        """Visualise set of hotspots using Lovis4u.

        Arguments:
            hotspot_ids (list): List of hotspot ids to be plotted.
            output_folder (str): Output folder to save pdf file.
            island_ids (None | list): List of island ids. In case it's specified only listed islands will be plotted.
            proteome_ids (None | list): List of proteome ids. In case it's specifiedd only listed proteome will
                be plotted.
            additional_annotation (dict): Additional LoVis4u feature annotation dict.
            keep_while_deduplication (list): List of island ids to be kept during deduplication.
            shortest_labels (bool| auto): Whether to put 1-based cds index only instead of CDS id for hypothetical
                proteins.

        Returns:
            None

        """
        try:
            if output_folder == "default":
                output_folder = os.path.join(self.prms.args["output_dir"], "lovis4u_hotspots")
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)
            locus_annotation_rows, feature_annotation_rows, mmseqs_results_rows, gff_files = [], [], [], []
            cds_table_rows = []
            already_added_groups = []
            hotspot_subset = self.hotspots.hotspots.loc[hotspot_ids].to_list()
            added_proteomes = []
            for hotspot in hotspot_subset:
                for h_island in hotspot.islands:
                    if island_ids:
                        if h_island.island_id not in island_ids:
                            continue
                    proteome = self.proteomes.proteomes.at[h_island.proteome]
                    if proteome_ids:
                        if proteome.proteome_id not in proteome_ids:
                            continue
                    if proteome.proteome_id in added_proteomes:
                        continue
                    added_proteomes.append(proteome.proteome_id)
                    proteome_annotation = self.proteomes.annotation.loc[h_island.proteome]
                    proteome_cdss = proteome.cdss.to_list()
                    locus_indexes = h_island.get_all_locus_indexes(proteome.cdss)
                    locus_groups = h_island.get_locus_groups(proteome.cdss)
                    if locus_groups not in already_added_groups:
                        already_added_groups.append(locus_groups)
                    else:
                        if h_island.island_id not in keep_while_deduplication:
                            continue
                    gff_files.append(proteome.gff_file)
                    start_coordinate = proteome_cdss[locus_indexes[0]].start
                    end_coordinate = proteome_cdss[locus_indexes[-1]].end
                    if compact_mode:
                        nf = 1
                        start_coordinate = proteome_cdss[h_island.indexes[0] - nf].start  #
                        end_coordinate = proteome_cdss[h_island.indexes[-1] + nf].end  #
                    if end_coordinate > start_coordinate:
                        sequence_coordinate = f"{start_coordinate}:{end_coordinate}:1"
                    else:
                        sequence_coordinate = f"{start_coordinate}:{proteome_annotation['length']}:1,1:{end_coordinate}:1"
                    locus_annotation_row = dict(sequence_id=h_island.proteome, coordinates=sequence_coordinate,
                                                circular=proteome.circular, group=hotspot.proteome_community)
                    if len(hotspot_ids) > 1:
                        locus_annotation_row["description"] = f"proteome community: {hotspot.proteome_community}"
                    locus_annotation_rows.append(locus_annotation_row)
                    for cds_ind, cds in enumerate(proteome_cdss):
                        if cds_ind in locus_indexes:
                            short_id = cds.cds_id.replace(proteome.proteome_id, "").strip().strip("_").strip("-")
                            if cds_ind not in h_island.indexes:
                                if cds.g_class == "conserved":
                                    group_type = "conserved"
                                    fcolour = "#8B9697"  # attention
                                else:
                                    group_type = "conserved"
                                    fcolour = "#D3D5D6"
                                scolour = "#000000"
                            else:
                                fcolour = "default"
                                scolour = "default"
                                group_type = "variable"
                                cds_table_row = dict(hotspot=hotspot.hotspot_id, sequence=h_island.proteome,
                                                     island_index=h_island.indexes.index(cds_ind),
                                                     group=cds.group,
                                                     coordinates=f"{cds.start}:{cds.end}:{cds.strand}",
                                                     cds_id=cds.cds_id, cds_type=cds.g_class,
                                                     name=cds.name)
                                if cds.hmmscan_results:
                                    cds_table_row["name"] = cds.hmmscan_results["target"]
                                    if "db_name" in cds.hmmscan_results.keys():
                                        cds_table_row["category"] = cds.hmmscan_results["db_name"]
                                    else:
                                        cds_table_row["category"] = cds.hmmscan_results["db"]
                                cds_table_rows.append(cds_table_row)
                            feature_annotation_row = dict(feature_id=cds.cds_id, group=cds.group, group_type=group_type,
                                                          fill_colour=fcolour, stroke_colour=scolour,
                                                          name=cds.name)
                            if feature_annotation_row["name"] == "hypothetical protein":
                                feature_annotation_row["name"] = ""
                            if cds.hmmscan_results:
                                feature_annotation_row["name"] = cds.hmmscan_results["target"]
                                if "db_name" in cds.hmmscan_results.keys():
                                    feature_annotation_row["category"] = cds.hmmscan_results["db_name"]
                                else:
                                    feature_annotation_row["category"] = cds.hmmscan_results["db"]
                            if additional_annotation:
                                if cds.cds_id in additional_annotation.keys():
                                    feature_annotation_row.update(additional_annotation[cds.cds_id])
                            if feature_annotation_row["name"]:
                                if not compact_mode:
                                    feature_annotation_row["name"] += f" ({short_id})"
                            else:
                                if shortest_labels == True or \
                                        (shortest_labels == "auto" and len(h_island.indexes) >= self.prms.args[
                                            "island_size_cutoff_to_show_index_only"]):
                                    feature_annotation_row["name"] = str(cds_ind + 1)
                                    if compact_mode:
                                        feature_annotation_row["name"] = ""
                                else:
                                    feature_annotation_row["name"] = str(short_id)
                            feature_annotation_rows.append(feature_annotation_row)
                        mmseqs_results_rows.append(dict(cluster=cds.group, protein_id=cds.cds_id))

            cds_tables_folder = os.path.join(self.prms.args["output_dir"], "lovis4u_hotspots_annotation")
            if not os.path.exists(cds_tables_folder):
                os.mkdir(cds_tables_folder)
            cds_table = pd.DataFrame(cds_table_rows)
            table_name = f"{'_'.join(hotspot_ids)}"
            if len(table_name) > 200:
                table_name = f"{table_name[:200]}..._{hotspot_ids[-1]}"
            cds_table.to_csv(os.path.join(cds_tables_folder, f"{table_name}.tsv"), sep="\t", index=False)

            locus_annotation_t = pd.DataFrame(locus_annotation_rows)
            feature_annotation_t = pd.DataFrame(feature_annotation_rows)
            temp_input_f = tempfile.NamedTemporaryFile()
            temp_input_l = tempfile.NamedTemporaryFile()
            locus_annotation_t.to_csv(temp_input_l.name, sep="\t", index=False)
            feature_annotation_t.to_csv(temp_input_f.name, sep="\t", index=False)

            l_parameters = lovis4u.Manager.Parameters()
            if compact_mode:
                l_parameters.load_config("A4p1")
            else:
                l_parameters.load_config(self.prms.args["lovis4u_hotspot_config_filename"])
            # l_parameters.load_config(self.prms.args["lovis4u_hotspot_config_filename"])
            l_parameters.args["cluster_all_proteins"] = False
            l_parameters.args["locus_label_style"] = "id"
            l_parameters.args["locus_label_position"] = "bottom"
            l_parameters.args["verbose"] = False
            l_parameters.args["draw_individual_x_axis"] = False

            l_parameters.args["draw_middle_line"] = True
            l_parameters.args["category_colours"] = self.prms.args["category_colours"]
            l_parameters.args["output_dir"] = os.path.join(self.prms.args["output_dir"], "lovis4u_tmp")

            if compact_mode:
                l_parameters.args["feature_group_types_to_show_label"] = []
                l_parameters.args["feature_group_types_to_show_label_on_first_occurrence"] = ["conserved", "variable"]
                l_parameters.args["draw_individual_x_axis"] = False

            loci = lovis4u.DataProcessing.Loci(parameters=l_parameters)
            loci.load_feature_annotation_file(temp_input_f.name)
            loci.load_locus_annotation_file(temp_input_l.name)

            mmseqs_results_t = pd.DataFrame(mmseqs_results_rows).set_index("protein_id")
            loci.load_loci_from_extended_gff(gff_files, ilund4u_mode=True)
            loci.cluster_sequences(mmseqs_results_t, one_cluster=True)
            loci.reorient_loci(ilund4u_mode=True)
            loci.set_feature_colours_based_on_groups()
            loci.set_category_colours()
            loci.define_labels_to_be_shown()

            loci.save_feature_annotation_table()
            loci.save_locus_annotation_table()

            canvas_manager = lovis4u.Manager.CanvasManager(l_parameters)
            canvas_manager.define_layout(loci)
            canvas_manager.add_loci_tracks(loci)
            canvas_manager.add_scale_line_track()
            canvas_manager.add_categories_colour_legend_track(loci)
            canvas_manager.add_homology_track()
            pdf_name = f"{'_'.join(hotspot_ids)}"
            if len(pdf_name) > 200:
                pdf_name = f"{pdf_name[:200]}..._{hotspot_ids[-1]}"
            canvas_manager.plot(f"{pdf_name}.pdf")
            os.system(f"mv {l_parameters.args['output_dir']}/{pdf_name}.pdf {output_folder}/")
            if keep_temp_data:
                if not os.path.exists(os.path.join(output_folder, "lovis4u_output")):
                    os.mkdir(os.path.join(output_folder, "lovis4u_output"))
                os.system(f"mv {l_parameters.args['output_dir']} "
                          f"{os.path.join(output_folder, 'lovis4u_output', pdf_name)}")
                os.mkdir(os.path.join(output_folder, "lovis4u_output", pdf_name, "gff_files"))
                for gff_file in gff_files:
                    os.system(f"cp '{gff_file}' {os.path.join(output_folder, 'lovis4u_output', pdf_name, 'gff_files')}/")
            else:
                shutil.rmtree(l_parameters.args["output_dir"])
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to plot set of hotspots using LoVis4u.") from error

    def plot_proteome_communities(self):
        """Run visualisation of proteome list for each proteome community.

        Returns:
            None

        """
        try:
            if self.prms.args["verbose"]:
                print(f"○ Visualisation of proteome communities with corresponding hotspots using lovis4u...",
                      file=sys.stdout)
                vis_output_folder = os.path.join(self.prms.args["output_dir"], "lovis4u_proteome_communities")
                if os.path.exists(vis_output_folder):
                    shutil.rmtree(vis_output_folder)
                os.mkdir(vis_output_folder)
                bar = progress.bar.FillingCirclesBar(" ", max=len(self.proteomes.communities.keys()),
                                                     suffix='%(index)d/%(max)d')
                for com_id, com_pr_ids in self.proteomes.communities.items():
                    bar.next()
                    self.plot_proteome_community(com_id, vis_output_folder)
                bar.finish()
        except Exception as error:
            raise ilund4u.manager.ilund4uError("Unable to plot proteome communities.") from error

    def plot_proteome_community(self, community: int, output_folder: str, mode: str = "hotspot",
                                proteome_ids: typing.Union[None, list] = None,
                                additional_annotation: typing.Union[None, dict] = None,
                                filename: typing.Union[None, str] = None):
        """Visualise proteome community using LoVis4u.

        Arguments:
            community (int): Id of proteome community to be plotted.
            output_folder (str): Output folder to save pdf file.
            mode (str): Mode of visualisation.
            proteome_ids: (None | list): List of proteome ids. In case it's specified only listed proteomes will
                be plotted.
            additional_annotation (dict): Additional LoVis4u feature annotation dict.
            filename (None | str): Pdf file name. If not specified id of community will be used.

        Returns:
            None

        """
        try:
            community_proteomes = self.proteomes.communities[community]
            if mode == "hotspot":
                hotspot_annotation_com = self.hotspots.annotation[
                    self.hotspots.annotation["proteome_community"] == community]
                num_of_hotspots = len(hotspot_annotation_com.index)
                colours_rgb = seaborn.color_palette("husl", num_of_hotspots, desat=1)
                colours = list(map(lambda x: matplotlib.colors.rgb2hex(x), colours_rgb))
                colours_dict = ({g: c for g, c in zip(list(hotspot_annotation_com.index.to_list()), colours)})
                com_hotspots = self.hotspots.hotspots.loc[hotspot_annotation_com.index]
                island_proteins_d = dict()
                for hotspot in com_hotspots.to_list():
                    for island in hotspot.islands:
                        proteome = self.proteomes.proteomes.at[island.proteome]
                        island_indexes = island.indexes
                        island_cds_ids = proteome.cdss.iloc[island_indexes].apply(lambda cds: cds.cds_id).to_list()
                        for ic_id in island_cds_ids:
                            island_proteins_d[ic_id] = hotspot.hotspot_id
            gff_files = []
            feature_annotation_rows = []
            mmseqs_results_rows = []
            n_of_added_proteomes = 0
            for proteome_id in community_proteomes:
                if proteome_ids:
                    if proteome_id not in proteome_ids:
                        continue
                n_of_added_proteomes += 1
                proteome = self.proteomes.proteomes.at[proteome_id]
                gff_files.append(proteome.gff_file)
                for cds_ind, cds in enumerate(proteome.cdss.to_list()):
                    if cds.g_class == "conserved":
                        group_type = "conserved"
                    else:
                        group_type = "variable"
                    if mode == "hotspot":
                        if cds.cds_id in island_proteins_d.keys():
                            fcolour = colours_dict[island_proteins_d[cds.cds_id]]
                        else:
                            if cds.g_class == "conserved":
                                fcolour = "#BDC6CA"
                            else:
                                fcolour = "#8C9295"
                    feature_annotation_row = dict(feature_id=cds.cds_id, group=cds.group, group_type=group_type)
                    if mode == "hotspot":
                        feature_annotation_row["show_label"] = 0
                        feature_annotation_row["stroke_colour"] = "#000000"
                        feature_annotation_row["fill_colour"] = fcolour
                    if cds.hmmscan_results and self.prms.args["show_hmmscan_hits_on_full_proteomes"]:
                        feature_annotation_row["name"] = cds.hmmscan_results["target"]
                        if "db_name" in cds.hmmscan_results.keys():
                            feature_annotation_row["category"] = cds.hmmscan_results["db_name"]
                        else:
                            feature_annotation_row["category"] = cds.hmmscan_results["db"].lower()
                    if additional_annotation:
                        if cds.cds_id in additional_annotation.keys():
                            feature_annotation_row.update(additional_annotation[cds.cds_id])
                    feature_annotation_rows.append(feature_annotation_row)
                    mmseqs_results_rows.append(dict(cluster=cds.group, protein_id=cds.cds_id))

            l_parameters = lovis4u.Manager.Parameters()
            l_parameters.load_config(self.prms.args["lovis4u_proteome_config_filename"])
            l_parameters.args["cluster_all_proteins"] = False
            if n_of_added_proteomes > 1:
                l_parameters.args["draw_individual_x_axis"] = False
            else:
                l_parameters.args["draw_individual_x_axis"] = True
            l_parameters.args["verbose"] = False
            l_parameters.args["locus_label_style"] = "id"
            if mode == "hotspot" and n_of_added_proteomes != 1:
                l_parameters.args["gff_CDS_category_source"] = "-"
            l_parameters.args["draw_middle_line"] = False
            l_parameters.args["category_colours"] = self.prms.args["category_colours"]
            l_parameters.args["output_dir"] = os.path.join(self.prms.args["output_dir"], "lovis4u_tmp")
            if os.path.exists(l_parameters.args["output_dir"]):
                shutil.rmtree(l_parameters.args["output_dir"])
            os.mkdir(l_parameters.args["output_dir"])
            loci = lovis4u.DataProcessing.Loci(parameters=l_parameters)
            feature_annotation_t = pd.DataFrame(feature_annotation_rows)
            temp_input_f = tempfile.NamedTemporaryFile()
            feature_annotation_t.to_csv(temp_input_f.name, sep="\t", index=False)
            loci.load_feature_annotation_file(temp_input_f.name)
            mmseqs_results_t = pd.DataFrame(mmseqs_results_rows).set_index("protein_id")
            loci.load_loci_from_extended_gff(gff_files, ilund4u_mode=True)
            if len(gff_files) <= self.prms.args["max_number_of_seqs_to_redefine_order"]:
                loci.cluster_sequences(mmseqs_results_t, one_cluster=True)
            loci.reorient_loci(ilund4u_mode=True)
            if mode == "regular" or n_of_added_proteomes == 1:
                loci.define_labels_to_be_shown()
            loci.set_feature_colours_based_on_groups()
            loci.set_category_colours()
            loci.save_feature_annotation_table()
            canvas_manager = lovis4u.Manager.CanvasManager(l_parameters)
            canvas_manager.define_layout(loci)
            canvas_manager.add_loci_tracks(loci)
            if n_of_added_proteomes > 1:
                canvas_manager.add_scale_line_track()
            canvas_manager.add_categories_colour_legend_track(loci)
            canvas_manager.add_homology_track()
            if not filename:
                filename = f"{community}.pdf"
            canvas_manager.plot(filename)
            os.system(f"mv {l_parameters.args['output_dir']}/{filename} {output_folder}/")
            shutil.rmtree(l_parameters.args["output_dir"])
            return None
        except Exception as error:
            raise ilund4u.manager.ilund4uError(f"Unable to plot proteome community {community}.") from error
