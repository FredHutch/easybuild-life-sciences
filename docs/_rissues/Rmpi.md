---
title: Rmpi Failures 
date: 2016-10-10
---

### MPI/Rmpi failures
We have seen multiple cases where an Rmpi job hangs for unknown reasons.  Often times all of the worker processes will be idle except for one, which will be running at 100% on one of the allocated cores.

We do not have a workaround for this at this time.  If you are using Rmpi you will either need to load R/3.2.2 as indicated above or re-factor your script to use a different paralelizing tool.

If you are having trouble with Rmpi based jobs using the latest R versions, we have two possible solutions.  First, ensure that the environment variable OPENBLAS_NUM_THREADS is set to 1.  OpenBLAS (a high-performance open-source basic linear algebra subprogram library) uses all available processors for some parallel operations.  This conflicts with multithreaded and MPI applications leading to problems with Rmpi scripts.  Setting that environment variable turns off OpenBLAS threading.

We have changed the module for the latest R versions to set this environment variable automatically.  However, if you think OpenBLAS threading may help your application, make sure to use the value of "SLURM_CPUS_ON_NODE" to set OPENBLAS_NUM_THREADS.

Second: we have discovered that some commonly used Rmpi functions do not appear to work correctly with the latest OpenMPI versions in the later R bundles. These are the "low level" Rmpi calls like rmpi.bcast.cmd() or rmpi.bcast() which should be replaced with calls to mpi.bcast.Robj2slave() and rmpi.bcast.Rfun2slave().
