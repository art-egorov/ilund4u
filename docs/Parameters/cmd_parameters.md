# Сommand-line parameters



## No mode

### Post-install steps

`--data`
:    Creates the 'ilund4u_data' folder in the current working directory.
     The folder contains adjustable configuration files used by iLund4u.

`--linux`
:    Replaces the mmseqs path in the pre-made config file from the MacOS
     version [default] to the Linux version.
`--mac`
:    Replaces the mmseqs path in the pre-made config file from the Linux
     version to the MacOS version.

`--get-hmms`
:    Download HMMs (hmmscan format) of defence, anti-defence, virulence,
     and AMR proteins from our server [data-sharing.atkinson-lab.com]

`--database <phages|plasmids>`
:    Download iLund4u database from our server [data-sharing.atkinson-lab.com].
     Available databases: "phages" and "plasmids" (see documentation for details)

### Miscellaneous arguments

`-h --help`
:    Show this help message and exit.

`-v, --version`
:    Show program version.

`--debug`
:    Provide detailed stack trace for debugging purposes.

`-q, --quiet`
:    Don't show progress messages.


**iLund4u can be used in three different modes:**

- `ilund4u hotspots [parameters]` - hotspot annotation.

- `ilund4u protein [parameters]` - protein search versus iLund4u database.

- `ilund4u proteome [parameters]` - proteome annotation and search versus iLund4u database.

To show mode-specific help messages use: `ilund4u --help [mode]`
(e.g. `ilund4u --help hotspots`)

---

## Hotspot annotation mode

### Mandatory arguments

`-gff <folder>`
:    Path to a folder containing extended gff files.
     Each gff file should contain a corresponding nucleotide sequence.
     (designed to handle pharokka/prokka produced annotation files).

### Optional arguments | data processing

`-ufid, --use-filename-as-id`
:    Use filename (wo extension) as contig id instead
     of the contig id written in a gff file.

`-mps, --min-proteome-size <int>`
:    Minimal number of proteins in a genome to be taken in analysis.
     [default: 15]

`-gct, --genome-circularity-table <path>`
:    Path to a table containing information for each genome about whether
     it is circular or not.
     Format: tab-separated; column names: id, circular. Values 1 or 0.
     For genomes which id is not listed, the default value will be used.
     [default: non circular; use `-cg/--circular-genomes` to change]

`--pclt, --protein-clusters-table <path>`
:    Path to a table with predefined protein clusters (to replace
    mmseqs clustering). Format: mmseqs-like cluster table
    (two columns, first - cluster_id; second - protein_id; without header)

`-psc, --proteome-sim-cutoff <float>`
:    Minimal fraction of shared homologous proteins between two
    proteomes to be connected in proteome network by edge.
    [default: 0.65]

`--pcot, --proteome-communities-table <path>`
:    Path to a table with predefined proteome communities (to replace network-based
    Leiden community detection algorithm). Format: mmseqs-like cluster table
    (two columns, first - community_id; second - proteome_id; without header)

`-spc, --single-proteome-community`
:    Set the same single community for all proteomes. Could be considered
     as pangenome mode when input consists of the same or highly similar species.

`-mpcs, --min-proteome-community-size <int>`
:   Minimal size of proteome community which then will be taken in
   hotspot search analysis [default: 10].

`-vpc, --variable-protein-cutoff <float>`
:    Upper limit of fraction of proteomes in which a protein homology
    group should be found within a proteome community to be called
    "variable". [default: 0.25]

`-cpc, --conserved-protein-cutoff <float>`
:    Lower limit of fraction of proteomes in which a protein homology
    group should be found within a proteome community to be called
    "conserved". [default: 0.75]

`-cg, --circular-genomes`
:    Consider each input genome as circular. That means that first
    and last proteins will be considered as neighbours.
    [default: non circular]

`-ncg, --non-circular-genomes`
:    Consider each input genome as non-circular. That means that first
    and last proteins won't be considered as neighbours.
    [default: circular]

`-hpc, --hotspot-presence-cutoff`
:    Fraction of proteomes within a proteome community in which clustered
    variable islands should be found to be considered as a hotspot.
    [default: 0.3]

`-rnf, --report-not-flanked`
:    Report in results hotspots that have flanked conserved genes only on one
    side (located on the end of non-circular sequences) [default: False]

### Optional arguments | others

`-o <dir name>`
:    Output dir name. This will be created if it does not exist.
	[default: ilund4u_{current_date}; e.g. ilund4u_2022_07_25-20_41]

`-o-db, --output-database <dirname>`
:    Output dir name for the database.
    If not specified then no database will be saved.

`-c <standard|<file.cfg>`
:    Path to a configuration file or name of a premade config file
    [default: standard]


### Miscellaneous arguments

`--debug`
:    Provide a detailed stack trace for debugging purposes.

`-q, --quiet`
:    Don't show progress messages.


---

## Protein annotation mode


### Mandatory arguments


`-fa <path>`
:    Path to a fasta file with query protein sequence.

`-db <path>`
:    Path to iLund4u database.

### Optional arguments | data processing

`-hsm, --homology-search-mode <group|proteins>`
:    Mode to define homologous proteins from the database.
    If "group" is selected, a query protein is assigned to the
    same homologous group as that of the best search hit
    protein. In "proteins" mode, all target proteins that
    pass the cutoffs are considered to be homologues.
    [default: group]

`-msqc, --mmseqs-query-cov <float>`
:    MMseqs search query coverage cutoff [default: 0.65]

`-mstc, --mmseqs-target-cov <float>`
:    MMseqs search query coverage cutoff [default: 0.65]

`-msf, --mmseqs-fident <float>`
:    MMseqs search fident (fraction of identical matches) cutoff
    [default: 0.3]

`-mse, --mmseqs-evalue`
 :   MMseqs search evalue cutoff [default: 1e-5]

`-fm, --fast-mmseqs`
:    MMseqs search only against a database of representative sequences.
     [default: False]

`-rnf, --report-not-flanked`
:    Report in results hotspots that have flanked conserved genes only on one
    side (located on the end of non-circular sequences) [default: False]

`-ql, --query-label <str>`
 :   Label for query protein homologues on LoVis4u visualisation.
    [default: database protein names]

### Optional arguments | others

`-o <name>`
 :   Output dir name. This will be created if it does not exist.
	[default: ilund4u_{current_date}; e.g. ilund4u_2022_07_25-20_41]

`-c <standard|<file.cfg>`
 :   Path to a configuration file or name of a premade config file
    [default: standard]


### Miscellaneous arguments

`--debug`
:    Provide a detailed stack trace for debugging purposes.

`-q, --quiet`
:    Don't show progress messages.

---

## Proteome annotation mode

### Mandatory arguments

`-gff <path>`
 :   Path to a query proteome gff file.
    (designed to handle pharokka/prokka produced annotation files).

`-db <path>`
 :   Path to iLund4u database.

### Optional arguments | data processing

`-ncg, --non-circular-genomes`
 :   Consider tje query genome as non-circular. This means that first
    and last proteins won't be considered to be neighbours.
    [default: circular]

**Note**: MMseqs arguments below are used for assigning homology
      group for each protein in your query proteome.

`-msqc, --mmseqs-query-cov <float>`
:   MMseqs search query coverage cutoff [default: 0.65]

`-mstc, --mmseqs-target-cov <float>`
 :   MMseqs search query coverage cutoff [default: 0.65]

`-msf, --mmseqs-fident <float>`
:    MMseqs search fident (fraction of identical matches) cutoff
    [default: 0.3]

`-mse, --mmseqs-evalue`
 :   MMseqs search evalue cutoff [default: 1e-5]

`-fm, --fast-mmseqs`
:    MMseqs search only against a database of representative sequences.
     [default: False]

`-rnf, --report-not-flanked`
:    Report in results hotspots that have flanked conserved genes only on one
    side (located on the end of non-circular sequences) [default: False]


### Optional arguments | others

`-o <name>`
 :   Output dir name. It will be created if it does not exist.
	[default: ilund4u_{current_date}; e.g. ilund4u_2022_07_25-20_41]

`-c <standard|<file.cfg>`
 :   Path to a configuration file or name of a premade config file
    [default: standard]


### Miscellaneous arguments

`--debug`
:    Provide a detailed stack trace for debugging purposes.

`--parsing-debug`
:    Provide detailed stack trace for debugging purposes
     for failed reading of gff files.

`-q, --quiet`
:    Don't show progress messages.


---

