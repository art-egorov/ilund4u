"""
This module provides some methods (e.g. colour transformation, data copying, etc.) used by the tool.
"""
import matplotlib.colors
import shutil
import sys
import os
import re

import pandas as pd
import pyhmmer.plan7
import pyhmmer.easel
import pyhmmer

import progress.bar
import requests
import tarfile

import ilund4u


def adjust_paths(to: str) -> None:
    """Change paths in the internal config files for linux or mac.

    Arguments:
        to (str): mac | linux

    Returns:
        None

    """
    internal_dir = os.path.join(os.path.dirname(__file__), "ilund4u_data")
    config_files = ["standard.cfg"]
    for config_file in config_files:
        config_file_path = os.path.join(internal_dir, config_file)
        with open(config_file_path, "r+") as config:
            if to == "linux":
                if not os.path.exists(os.path.join(internal_dir, "bin/mmseqs_linux")):
                    os.system(f"unzip -q -d {os.path.join(internal_dir, 'bin/')} "
                              f"{os.path.join(internal_dir, 'bin/mmseqs_linux.zip')}")
                config_txt = re.sub(r"mmseqs_mac/bin/mmseqs", "mmseqs_linux/bin/mmseqs", config.read())
            else:
                config_txt = re.sub(r"mmseqs_linux/bin/mmseqs", "mmseqs_mac/bin/mmseqs", config.read())
            config.seek(0)
            config.truncate()
            config.write(config_txt)
    print(f"⦿ mmseqs path was adjusted to {to}", file=sys.stdout)
    return None


def copy_package_data() -> None:
    """Copy the ilund4u package data folder to your current dir.

    Returns:
        None

    """
    try:
        users_dir = os.path.join(os.getcwd(), "ilund4u_data")
        internal_dir = os.path.join(os.path.dirname(__file__), "ilund4u_data")
        if os.path.exists(users_dir):
            print("Warning: ilund4u_data folder already exists. Remove it or change its name first before "
                  "updating with default.")
            return None
        shutil.copytree(internal_dir, users_dir, ignore=shutil.ignore_patterns("help*", ".*", "HMMs*", "bin"))
        print("⦿ ilund4u_data folder was copied to the current working directory", file=sys.stdout)
        return None
    except Exception as error:
        raise ilund4u.manager.ilund4uError(f"Unable to copy ilund4u folder in your working dir.") from error


def get_color(name: str, parameters: dict) -> str:
    """Get HEX color by its name

    Arguments:
        name (str): name of a color.
        parameters (dict): Parameters' object dict.

    Returns:
        str: HEX color.

    """
    hex_c = parameters.args["palette"][parameters.args[name]]
    return hex_c


def get_colour_rgba(name: str, parameters: dict) -> tuple:
    """Get rgba colour by its name

    Arguments:
        name (str): name of a colour.
        parameters (dict): Parameters' object dict.

    Returns:
        tuple: RGBA colour

    """
    return *matplotlib.colors.hex2color(get_color(name, parameters)), parameters.args[f"{name}_alpha"]


def update_path_extension(path: str, new_extension: str) -> str:
    """Get path basename and replace its extension

    Arguments:
        path (str): path to a file
        new_extension (str): new extension

    Returns:
        str: basename of a file with new extension.

    """
    updated_filename = f"{os.path.splitext(os.path.basename(path))[0]}.{new_extension}"
    return updated_filename


def run_pyhmmer(query_fasta: str, query_size: int, prms: ilund4u.manager.Parameters) -> pd.DataFrame:
    """Run pyhmmer hmmscan for a set of query proteins

    Arguments:
        query_fasta (str): Path to a query fasta file.
        query_size (int): Number of query proteins.
        prms (ilund4u.manager.Parameters): Parameters class object that holds all arguments.

    Returns:
        pd.DataFrame: Table with hmmscan search results.

    """
    with pyhmmer.easel.SequenceFile(query_fasta, digital=True) as seqs_file:
        query_proteins = seqs_file.read_block()
    num_of_query_proteins = query_size

    hmmscan_output_folder = os.path.join(prms.args["output_dir"], "hmmscan")
    if os.path.exists(hmmscan_output_folder):
        shutil.rmtree(hmmscan_output_folder)
    os.mkdir(hmmscan_output_folder)

    databases_cname = ["hmm_defence", "hmm_virulence", "hmm_anti_defence", "hmm_amr"]
    databases_names = ["DefenceFinder and CasFinder Databases", "Virulence Factor Database",
                       "Anti-Prokaryotic Immune Systems Database (dbAPIS)", "AMRFinderPlus Database"]
    databases_short_names = ["Defence", "Virulence", "Anti-defence", "AMR"]

    alignment_table_rows = []
    for db_ind, db_name in enumerate(databases_cname):
        db_alignment_table_rows = []
        db_shortname = databases_short_names[db_ind]
        db_path = prms.args[db_name]
        db_full_name = databases_names[db_ind]
        if not os.path.exists(db_path):
            print(f"  ⦿ Database {db_full_name} was not found.", file=sys.stdout)
            continue
        hmm_files = [fp for fp in os.listdir(db_path) if os.path.splitext(fp)[1].lower() == ".hmm" and fp[0] != "."]
        hmms = []
        for hmm_file in hmm_files:
            hmms.append(pyhmmer.plan7.HMMFile(os.path.join(db_path, hmm_file)).read())
        if prms.args["verbose"]:
            print(f"  ⦿ Running pyhmmer hmmscan versus {db_full_name}...", file=sys.stdout)
            bar = progress.bar.FillingCirclesBar("   ", max=num_of_query_proteins, suffix="%(index)d/%(max)d")
        for hits in pyhmmer.hmmscan(query_proteins, hmms, E=1e-3, cpus=0):
            if prms.args["verbose"]:
                bar.next()
            for hit in hits:
                if hit.included:
                    for domain in hit.domains.reported:
                        if domain.i_evalue < 1e-3:
                            alignment = domain.alignment
                            alignment_row = dict(query=alignment.target_name.decode(),
                                                 target_db=db_shortname, target=hit.name.decode(),
                                                 hit_evalue=hit.evalue, di_evalue=domain.i_evalue,
                                                 q_from=alignment.target_from, q_to=alignment.target_to,
                                                 qlen=alignment.target_length, t_from=alignment.hmm_from,
                                                 t_to=alignment.hmm_to, tlen=alignment.hmm_length)
                            alignment_row["q_cov"] = round((alignment_row["q_to"] - alignment_row["q_from"]) / \
                                                           alignment_row["qlen"], 2)
                            alignment_row["t_cov"] = round((alignment_row["t_to"] - alignment_row["t_from"]) / \
                                                           alignment_row["tlen"], 2)
                            if alignment_row["q_cov"] >= prms.args["hmmscan_query_coverage_cutoff"] and \
                                    alignment_row["t_cov"] >= prms.args["hmmscan_hmm_coverage_cutoff"]:
                                if hit.description:
                                    alignment_row["t_description"] = hit.description.decode()
                                else:
                                    alignment_row["t_description"] = hit.name.decode()
                                db_alignment_table_rows.append(alignment_row)
        if prms.args["verbose"]:
            bar.finish()
        alignment_table_rows += db_alignment_table_rows
        db_alignment_table = pd.DataFrame(db_alignment_table_rows)
        if not db_alignment_table.empty:
            db_alignment_table = db_alignment_table.sort_values(by="hit_evalue", ascending=True)
            db_alignment_table = db_alignment_table.drop_duplicates(subset="query", keep="first").set_index("query")
            db_alignment_table.to_csv(os.path.join(hmmscan_output_folder, f"{db_shortname}.tsv"), sep="\t",
                                      index_label="query")
        n_hits = len(db_alignment_table.index)
        if prms.args["verbose"]:
            print(f"    Number of hits: {n_hits}", file=sys.stdout)
    alignment_table = pd.DataFrame(alignment_table_rows)
    if not alignment_table.empty:
        alignment_table = alignment_table.sort_values(by="hit_evalue", ascending=True)
        alignment_table = alignment_table.drop_duplicates(subset="query", keep="first").set_index("query")

    return alignment_table


def download_file_with_progress(url: str, local_folder: str) -> None:
    """Function for downloading a particular file from a web server.

    Arguments:
        url (str): Link to the file.
        local_folder (str): Path to a folder where file will be saved.

    Returns:
        None

    """
    try:
        response = requests.head(url)
        file_size = int(response.headers.get('content-length', 0))
        # Extract the original file name from the URL
        file_name = os.path.basename(url)
        local_path = os.path.join(local_folder, file_name)
        # Stream the file download and show progress bar
        with requests.get(url, stream=True) as r, open(local_path, 'wb') as f:
            bar = progress.bar.FillingCirclesBar(" ", max=file_size // 8192, suffix='%(percent)d%%')
            downloaded_size = 0
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    downloaded_size += len(chunk)
                    f.write(chunk)
                    bar.next()
            bar.finish()
        # Verify that the file was fully downloaded
        if downloaded_size != file_size:
            raise ilund4u.manager.ilund4uError(f"Downloaded file size ({downloaded_size} bytes) does not match "
                                               f"expected size ({file_size} bytes).")
        print(f"⦿ File was saved to {local_path}")
        if file_name.endswith('.tar.gz'):
            with tarfile.open(local_path, 'r:gz') as tar:
                tar.extractall(path=local_folder)
            print(f"⦿ Folder was successfully unarchived")
            os.remove(local_path)
    except Exception as error:
        raise ilund4u.manager.ilund4uError(f"Unable to get file from the {url}.") from error


def get_HMM_models(parameters) -> None:
    """Download HMM models

    Returns:
        None

    """
    try:
        url =  parameters["hmm_models"]
        internal_dir = os.path.join(os.path.dirname(__file__), "ilund4u_data")
        if os.path.exists(os.path.join(internal_dir, "HMMs")):
            print(f"○ HMMs folder already exists and will be rewritten...", file=sys.stdout)
        # Add checking if it's already downloaded
        print(f"○ Downloading HMM models...\n"
              f"  Source: {url}", file=sys.stdout)
        download_file_with_progress(url, internal_dir)
        return None
    except Exception as error:
        raise ilund4u.manager.ilund4uError(f"Unable to download HMM models.") from error

def get_ilund4u_db(parameters, db, path = "./") -> None:
    """Download ilund4u database

    Returns:
        None

    """
    try:
        url = parameters[f"{db}_db"]
        # Add checking if it's already downloaded
        print(f"○ Downloading iLund4u {db} database...\n"
              f"  Source: {url}", file=sys.stdout)
        download_file_with_progress(url, path)
        return None
    except Exception as error:
        raise ilund4u.manager.ilund4uError(f"Unable to download {db} database.") from error