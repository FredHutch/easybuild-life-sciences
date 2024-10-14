#!/bin/bash

ml Python/3.10.8-GCCcore-12.2.0
ml SciPy-bundle/2023.02-gfbf-2022b
ml matplotlib/3.7.0-gfbf-2022b
ml tqdm/4.64.1-GCCcore-12.2.0
ml scikit-learn/1.2.1-gfbf-2022b
ml MACS2/2.2.9.1-foss-2022b
ml h5py/3.8.0-foss-2022b
ml plotly.py/5.13.1-GCCcore-12.2.0
ml lxml/4.9.2-GCCcore-12.2.0
ml networkx/3.0-gfbf-2022b
ml pyBigWig/0.3.22-foss-2022b
ml Seaborn/0.12.2-foss-2022b
ml scanpy/1.9.3-foss-2022b
ml jupyter-server/2.7.0-GCCcore-12.2.0

export PYTHONPATH=/app/software/scenicplus/1.0.0-foss-2022b/lib/python3.10/site-packages:$PYTHONPATH

