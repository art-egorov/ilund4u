"""
This module provides managing classes and methods for the tool.
"""
import collections
import argparse
import configs
import time
import math
import sys
import os
import re
import json
import shutil

import reportlab.pdfbase.pdfmetrics as pdfmetrics
from reportlab.lib.units import cm, mm
import reportlab.pdfbase.ttfonts
import reportlab.pdfgen.canvas
import reportlab.rl_config
import reportlab.pdfgen

import ilund4u


class ilund4uError(Exception):
    """A class for exceptions parsing inherited from the Exception class.

    """
    pass


class Parameters:
    """A Parameters object holds and parse command line and config arguments.

    A Parameters object have to be created in each script since it's used almost by each
        class of the tool as a mandatory argument.

    Attributes:
        args (dict): dictionary that holds all arguments.
        cmd_arguments (dict): dictionary wich command-line arguments.

    """

    def __init__(self):
        """Parameters class constructor.

        """
        self.args = dict(debug=True)
        self.cmd_arguments = dict()

    def parse_cmd_arguments(self) -> None:
        """Parse command-line arguments.

        Returns:
            None

        """
        parser = argparse.ArgumentParser(prog="ilund4u", add_help=False)
        parser.add_argument("-data", "--data", dest="ilund4u_data", action="store_true")
        parser.add_argument("-get-hmms", "--get-hmms", dest="get_hmms", action="store_true")
        parser.add_argument("-database", "--database", dest="get_database", default=None, type=str,
                            choices=["phages", "plasmids"])

        parser.add_argument("-linux", "--linux", dest="linux", action="store_true", default=None)
        parser.add_argument("-mac", "--mac", dest="mac", action="store_true", default=None)
        parser.add_argument("-h", "--help", dest="help", action="store_true")
        parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.8")

        subparsers = parser.add_subparsers(dest="mode")

        parser_hm = subparsers.add_parser("hotspots")
        parser_hm.add_argument("-gff", "--gff", dest="gff", type=str, default=None)
        parser_hm.add_argument("-ufid", "--use-filename-as-id", dest="use_filename_as_contig_id", action="store_true",
                               default=None)
        parser_hm.add_argument("-mps", "--min-proteome-size", dest="min_proteome_size", type=int, default=None)
        parser_hm.add_argument("-gct", "--genome-circularity-table", dest="genome_annotation", type=str, default=None)
        parser_hm.add_argument("-psc", "--proteome-sim-cutoff", dest="proteome_similarity_cutoff", type=float,
                               default=None)
        parser_hm.add_argument("-mpcs", "--min-proteome-community-size", dest="min_proteome_community_size", type=int,
                               default=None)
        parser_hm.add_argument("-vpc", "--variable-protein-cutoff", dest="variable_protein_cluster_cutoff", type=float,
                               default=None)
        parser_hm.add_argument("-cpc", "--conserved-protein-cutoff", dest="conserved_protein_cluster_cutoff",
                               type=float, default=None)
        parser_hm.add_argument("-cg", "--circular-genomes", dest="circular_genomes", default=None,
                               action="store_true")
        parser_hm.add_argument("-ncg", "--non-circular-genomes", dest="circular_genomes", default=None,
                               action="store_false")
        parser_hm.add_argument("-hpc", "--hotspot-presence-cutoff", dest="hotspot_presence_cutoff",
                               type=float, default=None)
        parser_hm.add_argument("-rnf", "--report-not-flanked", dest="report_not_flanked", action="store_true",
                               default=None)
        parser_hm.add_argument("-o-db", "--output-database", dest="output_database", type=str, default=None)
        parser_hm.add_argument("-o", dest="output_dir", type=str, default=None)
        parser_hm.add_argument("-c", dest="config_file", type=str, default="standard")
        parser_hm.add_argument("-q", "--quiet", dest="verbose", default=True, action="store_false")
        parser_hm.add_argument("--debug", "-debug", dest="debug", action="store_true")
        parser_hm.add_argument("--parsing-debug", "-parsing-debug", dest="parsing_debug", action="store_true")

        parser_ps = subparsers.add_parser("protein")
        parser_ps.add_argument("-fa", "--fa", dest="fa", type=str, default=None)
        parser_ps.add_argument("-db", "--database", dest="database", type=str, default=None)
        parser_ps.add_argument("-ql", "--query-label", dest="query_label", type=str, default=None)
        parser_ps.add_argument("-hsm", "--homology-search-mode", dest="protein_search_target_mode", type=str,
                               choices=["group", "proteins"], default="group")
        parser_ps.add_argument("-msqc", "--mmseqs-query-cov", dest="mmseqs_search_qcov", type=float, default=None)
        parser_ps.add_argument("-mstc", "--mmseqs-target-cov", dest="mmseqs_search_tcov", type=float, default=None)
        parser_ps.add_argument("-msf", "--mmseqs-fident", dest="mmseqs_search_fident", type=float, default=None)
        parser_ps.add_argument("-mse", "--mmseqs-evalue", dest="mmseqs_search_evalue", type=float, default=None)
        parser_ps.add_argument("-fm", "--fast-mmseqs", dest="fast_mmseqs_search_mode", action="store_true",
                               default=None)
        parser_ps.add_argument("-rnf", "--report-not-flanked", dest="report_not_flanked", action="store_true",
                               default=None)
        parser_ps.add_argument("-c", dest="config_file", type=str, default="standard")
        parser_ps.add_argument("-o", dest="output_dir", type=str, default=None)
        parser_ps.add_argument("-q", "--quiet", dest="verbose", default=True, action="store_false")
        parser_ps.add_argument("--debug", "-debug", dest="debug", action="store_true")

        parser_pa = subparsers.add_parser("proteome")
        parser_pa.add_argument("-gff", "--gff", dest="gff", type=str, default=None)
        parser_pa.add_argument("-db", "--database", dest="database", type=str, default=None)
        parser_pa.add_argument("-ncg", "--non-circular-genomes", dest="circular_genomes", default=True,
                               action="store_false")
        parser_pa.add_argument("-msqc", "--mmseqs-query-cov", dest="mmseqs_search_qcov", type=float, default=None)
        parser_pa.add_argument("-mstc", "--mmseqs-target-cov", dest="mmseqs_search_tcov", type=float, default=None)
        parser_pa.add_argument("-msf", "--mmseqs-fident", dest="mmseqs_search_fident", type=float, default=None)
        parser_pa.add_argument("-mse", "--mmseqs-evalue", dest="mmseqs_search_evalue", type=float, default=None)
        parser_pa.add_argument("-fm", "--fast-mmseqs", dest="fast_mmseqs_search_mode", action="store_true",
                               default=None)
        parser_pa.add_argument("-rnf", "--report-not-flanked", dest="report_not_flanked", action="store_true",
                               default=None)
        parser_pa.add_argument("-c", dest="config_file", type=str, default="standard")
        parser_pa.add_argument("-o", dest="output_dir", type=str, default=None)
        parser_pa.add_argument("-q", "--quiet", dest="verbose", default=True, action="store_false")
        parser_pa.add_argument("--debug", "-debug", dest="debug", action="store_true")

        args = vars(parser.parse_args())
        if len(sys.argv[1:]) == 0:
            args["help"] = True
        if args["ilund4u_data"]:
            ilund4u.methods.copy_package_data()
            sys.exit()
        if args["linux"]:
            ilund4u.methods.adjust_paths("linux")
            sys.exit()
        if args["mac"]:
            ilund4u.methods.adjust_paths("mac")
            sys.exit()
        if args["get_hmms"]:
            self.load_config()
            ilund4u.methods.get_HMM_models(self.args)
            sys.exit()
        if args["get_database"]:
            self.load_config()
            ilund4u.methods.get_ilund4u_db(self.args, args["get_database"])
            sys.exit()
        if args["help"]:
            if not args["mode"]:
                help_message_path = os.path.join(os.path.dirname(__file__), "ilund4u_data", "help_main.txt")
            else:
                help_message_path = os.path.join(os.path.dirname(__file__), "ilund4u_data", f"help_{args['mode']}.txt")
            with open(help_message_path, "r") as help_message:
                print(help_message.read(), file=sys.stdout)
                sys.exit()
        if args["mode"] == "hotspots":
            if not args["gff"]:
                raise ilund4uError("-gff argument is required for hotspots mode.")
        elif args["mode"] == "protein":
            if not args["fa"] or not args["database"]:
                raise ilund4uError("Both -fa/--fa and -db/--database arguments are required for protein mode.")
        elif args["mode"] == "proteome":
            if not args["gff"] or not args["database"]:
                raise ilund4uError("Both -gff/--gff and -db/--database arguments are required for protein mode.")

        args_to_keep = ["gff", "output_database", "query_label", "genome_annotation"]
        filtered_args = {k: v for k, v in args.items() if v is not None or k in args_to_keep}
        self.cmd_arguments = filtered_args
        return None

    def load_config(self, path: str = "standard") -> None:
        """Load configuration file.

        Arguments
            path (str): path to a config file or name (only standard available at this moment).

        Returns:
            None

        """
        try:
            if path == "standard":
                path = os.path.join(os.path.dirname(__file__), "ilund4u_data", "standard.cfg")
            config = configs.load(path).get_config()
            internal_dir = os.path.dirname(__file__)
            for key in config["root"].keys():
                if type(config["root"][key]) is str and "{internal}" in config["root"][key]:
                    config["root"][key] = config["root"][key].replace("{internal}",
                                                                      os.path.join(internal_dir, "ilund4u_data"))
            config["root"]["output_dir"] = config["root"]["output_dir"].replace("{current_date}",
                                                                                time.strftime("%Y_%m_%d-%H_%M"))
            keys_to_transform_to_list = []
            for ktl in keys_to_transform_to_list:
                if isinstance(config["root"][ktl], str):
                    if config["root"][ktl] != "None":
                        config["root"][ktl] = [config["root"][ktl]]
                    else:
                        config["root"][ktl] = []
            self.args.update(config["root"])
            if self.cmd_arguments:
                self.args.update(self.cmd_arguments)
            return None
        except Exception as error:
            raise ilund4uError("Unable to parse the specified config file. Please check your config file "
                               "or provided name.") from error
