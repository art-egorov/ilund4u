# FAQ


 __*Q:* Something weird is happening, e.g. program is crashing with my data. What should I do?__  
 *A:* Please add `--debug` to the command line (or `--parsing-debug` if some files were not read correctly) and post the details in [Issues](https://github.com/art-egorov/ilund4u/issues).


 __*Q:* How to get representative protein ids of flanking genes for a hotspot?__  
 *A:* Each hotspot includes an attribute called *conserved signature*, which contains the protein IDs of representative sequences for protein families that are highly conserved as flanking genes (found as flanking genes in over 75% of hotspot islands). Starting with iLund4u version 0.0.8, this information is available in the *conserved_signature* column of the *hotspot_annotation.tsv* file.

 __*Q:* How to use pre-defined protein clusters instead of mmseqs clustering step?__  
 *A:* Starting with 0.0.10 version iLund4u has `--pclt, --protein-clusters-table <path>` parameter that takes as input path to a table with predefined protein clusters. Format: mmseqs-like cluster table (two columns, first - cluster_id; second - protein_id; without header).


 __*Q:* How to use pre-defined proteome (sequence/contigs) communities instead of network-based Leiden community detection algorithm?__  
 *A:* Starting with 0.0.10 version iLund4u has `--pcot, --proteome-communities-table <path>` parameter that takes as input path to a table with predefined proteome communities. Format: mmseqs-like cluster table (two columns, first - community_id; second - proteome_id; without header). In addition, it's possible to force iLund4u to consider all contigs as members of one community (pangenome-like mode). For that, simply apply `-spc, --single-proteome-community` parameter.
