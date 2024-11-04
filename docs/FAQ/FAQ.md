# FAQ


 __*Q:* Something weird is happening, e.g. program is crashing with my data. What should I do?__  
 *A:* Please add `--debug` to the command line (or `--parsing-debug` if some files were not read correctly) and post the details in [Issues](https://github.com/art-egorov/ilund4u/issues).


 __*Q:* How to get representative protein ids of flanking genes for a hotspot__  
 *A:* Each hotspot includes an attribute called *conserved signature*, which contains the protein IDs of representative sequences for protein families that are highly conserved as flanking genes (found as flanking genes in over 75% of hotspot islands). Starting with iLund4u version 0.0.8, this information is available in the *conserved_signature* column of the *hotspot_annotation.tsv* file.

