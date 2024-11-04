"""
This module provides data managing classes and methods for the tool.
"""
import shutil
import json
import time
import sys
import os

import ilund4u


class DatabaseManager:
    """Manager for loading and building iLund4u database.

    Attributes:
        prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

    """

    def __init__(self, parameters: ilund4u.manager.Parameters):
        """DatabaseManager class constructor.

        Arguments:
            parameters (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

        """
        self.prms = parameters

    def build_database(self, proteomes: ilund4u.data_processing.Proteomes, hotspots: ilund4u.data_processing.Hotspots,
                       db_path: str) -> None:
        """Write database.

        Arguments:
            proteomes (ilund4u.data_processing.Proteomes): Proteomes object.
            hotspots (ilund4u.data_processing.Hotspots): Hotspots object.
            db_path (str): Path to the database folder.

        Returns:

        """
        if os.path.exists(db_path):
            if self.prms.args["verbose"]:
                print("○ Warning: database folder will be rewritten.")
            shutil.rmtree(db_path)
        os.mkdir(db_path)
        if self.prms.args["verbose"]:
            print(f"○ Database building...", file=sys.stdout)
        proteomes.save_as_db(db_path)
        hotspots.save_as_db(db_path)
        database_info_txt = f"Date and time of building: {time.strftime('%Y.%m.%d-%H:%M')}\n" \
                            f"iLund4u version: {self.prms.args['version']}"
        with open(os.path.join(db_path, "db_info.txt"), "w") as db_info:
            db_info.write(database_info_txt)
        with open(os.path.join(db_path, "parameters.json"), "w") as parameters:
            json.dump(self.prms.args, parameters)
        if self.prms.args["verbose"]:
            print(f"  ⦿ Database was successfully saved to {db_path}", file=sys.stdout)
        return None

    def load_database(self, db_path: str) -> ilund4u.data_processing.Database:
        """Load database from its folder path and create a Database class object.

        Arguments:
            db_path (str): Path to the pre-built database folder.

        Returns:
            ilund4u.data_processing.Database: Database class object.

        """
        if self.prms.args["verbose"]:
            print(f"○ Loading database from {db_path}...", file=sys.stdout)
        proteomes = ilund4u.data_processing.Proteomes.db_init(db_path, self.prms)
        hotspots = ilund4u.data_processing.Hotspots.db_init(db_path, proteomes, self.prms)
        db_paths = dict(db_path=db_path, rep_fasta=os.path.join(db_path, "representative_seqs.fa"),
                        proteins_db=os.path.join(db_path, "mmseqs_db", "all_proteins"))
        if os.path.exists(os.path.join(db_path, "protein_group_accumulated_statistics.tsv")):
            db_paths["protein_group_stat"] = os.path.join(db_path, "protein_group_accumulated_statistics.tsv")

        database = ilund4u.data_processing.Database(proteomes, hotspots, db_paths, self.prms)
        if self.prms.args["verbose"]:
            print(f"⦿ The {db_path} database was successfully loaded", file=sys.stdout)
        return database
