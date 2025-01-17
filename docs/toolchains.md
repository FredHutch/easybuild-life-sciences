---
title: "Toolchains"
layout: single
permalink: /toolchains/
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

Most software packages are built from a compiler toolchain. A toolchain consists of one or more compilers, 
together with some libraries for using an MPI stack and commonly used math libraries. Well-known packages
include BLAS/LAPACK APIs for linear algebra routines.

EasyBuild supports several different toolchains. The Hutch uses the foss toolchain, which
is an acronym for Free Open Source Software. Foss is based on the GNU compiler tools.
Toolchains are curated twice a year and named for the year, and labeled a or b, i.e., foss-2019b.
When loading multiple modules for a workflow, always load modules from the same toolchain. Do not mix and
match modules from different toolchains.

## FOSS Toolchain Versions

The foss toolchain consists entirely of open source software. The FOSS name is derived from Free Open Sorce Software.

Each toolchain uses a single version of Python. Multible versions of libraries can exist within a single toolchain, like 
`matplotlib` but the base Python remains the same for all libraries within a toolchain.


| Toolchain | GCC Version | Python | CUDA  |Notes |
| ----------|-------------|--------|-------|------|
| [foss-2024a]({{ site.baseurl }}/toolchains/foss-2024a/) | GCC 13.3.0 | 3.12.3 | 12.6.2 |November 2024 Harmony |
| [foss-2023b]({{ site.baseurl }}/toolchains/foss-2023b/) | GCC 13.2.0 | 3.11.5 | 12.4.0 |March 19, 2024 |
| [foss-2023a]({{ site.baseurl }}/toolchains/foss-2023a/) | GCC 12.3.0 | 3.11.3 | 12.1.1 | Begining Feb 2024 |
| [foss-2022b]({{ site.baseurl }}/toolchains/foss-2022b/) | GCC 12.2.0 | 3.10.8 | 11.7.0 |Begining Feb 2023 |
| [foss-2022a]({{ site.baseurl }}/toolchains/foss-2021b/) | GCC 11.3.0 | 3.10.4 | | Limited use |
| [foss-2021b]({{ site.baseurl }}/toolchains/foss-2021b/) | GCC 11.2.0 | 3.9.6  | | Begining Nov 2021 |
| [foss-2021a]({{ site.baseurl }}/toolchains/foss-2021a/) | GCC 10.3.0 | 3.9.5  | | Limited use |
| [foss-2020b]({{ site.baseurl }}/toolchains/foss-2020b/) | GCC 10.2.0 | 3.8.6  | | Begining Nov 2020 |
| foss-2020a                                              | GCC 9.3.0  | 3.8.2  | | Mostly skipped |
| [foss-2019b]({{ site.baseurl }}/toolchains/foss-2019b/) | GCC 8.3.0  | 3.7.4  | | Primay tool chain for June 2019 cluster update |
| foss-2019a | GCC 8.2.0 |  | skipped |
| [foss-2018b]({{ site.baseurl }}/toolchains/foss-2018b/) | GCC 7.3.0 | | Many modules were built for the new cluster, but were outdated by 2019 |
| [foss-2016b]({{ site.baseurl }}/toolchains/foss-2016b/) | GCC 5.4.0 | | In use from 2016 too 2019. Retired with the Ubuntu 14.04 cluster |

## Toolchain diagram

To be more helpful in understanding the differences between these families, here is a diagram that explains what is added in
each additional layer.

Note: because there have been a few changes in toolchains, there are notes below the diagram
that explain the differences between the generations going back to the `2020b` version of the `foss` and `intel` toolchains.

<!-- https://github.com/easybuilders/easybuild-docs/blob/03891cbe6404a7fa237f289c99a660cfac5d7a73/docs/common-toolchains.md?plain=1#L9 -->
### Newest generations (`2023a` and later):


<!--

Mermaid diagrams will not render on GitHub Pages sites.
So I took a screenshot of the diagram as rendered on github.com
and display it here. If you need to update the diagram, uncomment 
the mermaid code below and update the diagram. Then take a screenshot
and save it as docs/images/toolchain-diagram.png.

Note that the Mermaid code has several "hyphen-hyphen greater than" in
it and that breaks the HTML comment. So I changed the hyphens to tildes,
if you modify the diagram, please change them back to hyphens.

-->

![Toolchain Diagram](images/toolchain-diagram.png)

<!--
```mermaid
graph LR
  A[GCCCore] ~~> |binutils| B[GCC];
  A ~~> |binutils| C[intel-compilers];
  B ~~> |OpenMPI| E[gompi];
  C ~~> |impi| F[iimpi];
  B ~~> |FlexiBLAS + FFTW + ScaLAPACK| D[gfbf];
  D ~~> |OpenMPI| G[foss];
  E ~~> |FlexiBLAS + FFTW + ScaLAPACK| G[foss];
  F ~~> |imkl| Z[intel];
  C ~~> |imkl| H[iimkl];
  H ~~> |impi| Z[intel];
```

-->

Note: following notes apply for the generations listed and those older than it:

- `2023a` - `iimkl` not present yet
- `2021b` - `gfbf` not present yet
- `2020b` - `foss` uses OpenBLAS instead of FlexiBLAS, `iccifort` is used instead of `intel-compilers`

