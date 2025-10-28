# Version log

* **Ver 0.1.4*** - 28 Octover
	- Incorrect parsing of pre-defined proteome community was fixed

* Ver 0.1.3 - 21 Octover
	- Hotspot visualisation but fix in protein search mode

* Ver 0.1.1 and 0.1.2 - 10 October
	- Update of the protein search mode. Now it can return results for multiple protein families that are homologous to the query sequence.

* Ver 0.1.0 - 17 September
	- Udpate of the default parameters and databases (the main - mmseqs clustering mmseqs_min_seq_id parameter updated from 0.3 to 0.5).

* Ver 0.0.11 - 26 January
	- Database reading was updated to cover pkl file parsing.

* **Ver 0.0.10** - 27 December
	- Parameters to use predefined protein and proteome clusters were introduced.
	- Log messages for protein mode were adjusted.
	- MMSeqs2 `--max-seqs` clustering parameter is now auto-adjusted for large input. 

* Ver 0.0.9 - 11 December
	- Now for dropped sequences information about duplication source is reported.
	- Config file parameter for mmseqs *--max-seqs* option is added.

* Ver 0.0.8.1 - 5 November
	- Bug with LoVis4u visualisation when `-ufid` parameter is used was fixed. 

* Ver 0.0.8 - 4 November
	- PADLOC hmm models were added.
	- Fast mode search for protein and proteome annotation mode is introduced.
	- Updated statistics output for all modes.

* [Preprint](https://doi.org/10.1101/2024.10.15.618418) released - 16 October,  ✨

* Ver 0.0.7.1 - 15 October 2024
	- `-use-filename-as-id` and `--parsing-deubg` parameters were added

* Ver 0.0.7 - 15 October 2024
	- Minor bug fixes
	- Visualisation adjusment

* Ver 0.0.6 - 30 September 2024 
	- Hotspots that are located on the end of non-circular chromosomes are now treated separately. 
	- `-rnf, --report-not-flanked` parameter is added.

* Ver 0.0.5.1 - 5 August 2024 
	- Long table file name bug is fixed (adjusment to the 0.0.5 v. update).

* Ver 0.0.5 - 25 July 2024 
	- Long pdf file name bug is fixed.

* Ver 0.0.4 - 27 June 2024 
	- Minor fixes and parameter name changes.

* Ver 0.0.3 - 27 June 2024 
	- Minor fixes.

* Ver 0.0.2 - 13 June 2024 
	- `-ncg`/`--non-circular-genomes` cmd parameter is added.
	- Minor fixes.

* Ver 0.0.1 - 13 June 2024 - first public release. 
