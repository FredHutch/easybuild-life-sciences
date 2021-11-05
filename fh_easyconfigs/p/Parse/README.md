build notes Nov, 4, 2021

Company is call ParseBioSciences
Python package name is splitpipe
Code is installed as split_seq-0.9.3p.dis-info and split_seq-0.9.3p-py3.8.egg
finished binary is named split-pipe
easybuild sanity check was fixed by using splitpipe as the package name.
but use split-pipe -h for sanity check command


louvain 0.7.0 was failing to install. setup.py was not seeing the package. It is trying to build igraph from source.
 - igraph lib locations are hardcoded into the setup.py
 - fhPython 3.8.6 has louvain 0.7.0

>>> import louvain
<stdin>:1: DeprecationWarning: This package has been superseded by the `leidenalg` package and will no longer be maintained. Please upgrade to the `leidenalg` package.
>>> louvain.__version__
'0.7.0'

Parse requires leidenalg, why does it also need Louvain?

remove Louvain and builds without issue.  
Louvain built with fhPython, How did I get that to work?

setup py is downlowning louvain, 0.7.0

Downloading https://files.pythonhosted.org/packages/e8/9d/dd6f02523ffbd12a2c396b7b82d91e60f6574d12648a190fe7756ce8bb66/louvain-0.7.0-cp38-cp38-manylinux2010_x86_64.whl

# test what happens when attempting to install louvain. The error message must break pip?
Requirement already satisfied: louvain in /app/software/Python/3.8.6-GCCcore-10.2.0/lib/python3.8/site-packages (0.7.0)
Requirement already satisfied: python-igraph>=0.8.0 in /app/software/python-igraph/0.9.0-foss-2020b/lib/python3.8/site-packages (from louvain) (0.9.0)
Requirement already satisfied: texttable>=1.6.2 in /app/software/python-igraph/0.9.0-foss-2020b/lib/python3.8/site-packages (from python-igraph>=0.8.0->louvain) (1.6.3)

