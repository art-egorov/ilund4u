import setuptools
import os

with open("docs/pypi.md", "r") as fh:
    long_description = fh.read()

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        if "ilund4u/ilund4u_data/HMMs" not in path:
            for filename in filenames:
                paths.append(os.path.join('..', path, filename))
    return paths


extra_files = package_files("ilund4u/ilund4u_data")
extra_files.append("../docs/pypi.md")

setuptools.setup(name="ilund4u",
                 version="0.0.8",
                 python_requires='>=3.8',
                 description="description",
                 url="https://art-egorov.github.io/ilund4u/",
                 author="Artyom Egorov",
                 author_email="artem.egorov@med.lu.se",
                 license="WTFPL",
                 packages=["ilund4u"],
                 package_data={"ilund4u": extra_files},
                 install_requires=["biopython", "requests", "configs", "argparse", "pandas", "bcbio-gff", "matplotlib",
                                   "seaborn", "scipy", "msa4u", "lovis4u", "progress", "leidenalg", "igraph", "pyhmmer"],
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 scripts=["bin/ilund4u"],
                 zip_safe=False)
