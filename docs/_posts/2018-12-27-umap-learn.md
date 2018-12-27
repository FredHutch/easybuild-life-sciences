---
layout: post
title: UMAP 
catagory:  
homepage: https://umap-learn.readthedocs.io/en/latest/ 
---
umap is installed as Python module.
Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction. The algorithm is founded on three assumptions about the data
 - The data is uniformly distributed on Riemannian manifold;
 - The Riemannian metric is locally constant (or can be approximated as such);
 - The manifold is locally connected.
From these assumptions it is possible to model the manifold with a fuzzy topological structure. The embedding is found by searching for a low dimensional projection of the data that has the closest possible equivalent fuzzy topological structure. 
```
module load Python/3.6.7-foss-2016b-fh2
```
