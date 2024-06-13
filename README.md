
<img  src="docs/img/ilund4u_logo.png" width="300"/>


## Description

**iLund4u** is a bioinformatics tool for search and annotation of hotspots in a large set of proteomes. 

**Supported input**: gff3 with sequence (prokka/pharokka produced)      
**Programming language:** Python3   
**OS:** MacOS, Linux  
**Python dependencies:** biopython, bcbio-gff, scipy, configs, argparse, pandas, matplotlib, seaborn, progess, leidanalg, igraph. pyhmmer, msa4u, lovis4u    
**Python version:** >= 3.8  
**OS-level dependencies:** MMseqs2 (included in the package)  
**License:** WTFPL  
**Version:** 0.0.1 (June 2024)

**Detailed documentation with user guide is available at [iLund4u Homepage](https://art-egorov.github.io/iLund4u/)**

<img  src="docs/img/ilundu4_pipeline.png" width="100%"/>


## Installation 🛠️

- The most stable release of ilund4u can be installed directly from pypi:

```
python3 -m pip install ilund4u
```

- The development version is available at github :

```
git clone https://github.com/art-egorov/ilund4u.git
cd ilund4u
python3 -m pip install --upgrade pip
python3 -m pip install setuptools wheel
python3 setup.py sdist
python3 -m pip install -e .
```

**!** If you're a linux user, run `ilund4u --linux` post-install command once to update paths in the premade config files that set by default for MacOS users.


## Reference 📃

If you find iLund4u useful, please cite:  
Artyom. A. Egorov, Gemma C. Atkinson, **iLund4u: ---**, *---*

## Contact 📇

Please contact us by e-mail _artem**dot**egorov**AT**med**dot**lu**dot**se_ or use [Issues](https://github.com/art-egorov/ilund4u/issues?q=) to report any technical problems.  
You can also use [Discussions section](https://github.com/art-egorov/ilund4u/discussions) for sharing your ideas or feature requests! 

## Authors 👨🏻‍💻

iLund4u is developed by Artyom Egorov at [the Atkinson Lab](https://atkinson-lab.com), Department of Experimental Medical Science, Lund University, Sweden. We are open for suggestions to extend and improve iLund4u functionality. Please don't hesitate to share your ideas or feature requests.
