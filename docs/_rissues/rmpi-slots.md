---
title: "Not enough slots" Starting Rmpi Cluster
date: 2016-01-01
---

With the recent updates of the R, Rmpi, and OpenMPI module bundles you may now see errors similar to this when running your R script:

```
There are not enough slots available in the system to satisfy the 3 slots 
that were requested by the application:
  /app/easybuild/software/R/3.3.3-foss-2016b-fh2/lib/R/library/Rmpi/Rslaves.sh

Either request fewer slots for your application, or make more slots available
for use.
```

This is because of a change in the interaction between Rmpi and OpenMPI.  It's now necessary to reduce the number of slaves spawned to account for the master process running the MPI session.  Typically your Rmpi script will contain something like this to spawn the slaves on the allocated machines:

```
library("Rmpi")
mpi.spawn.Rslaves()
...
```

The default for this call is to spawn one slave in every assigned slot which incorrectly assigns a slave into the master slot.  It is now necessary to reduce the number of slaves spawned by one:

```
library("Rmpi")
n=mpi.universe.size() - 1
mpi.spawn.Rslaves(nslaves=n)
...
```

this could, of course, be combined into a single line.  When successful, a simple "hello, world" like this:

```
#!/usr/bin/env Rscript

library('Rmpi')

# In this version we reduce the number of slaves started to account for the
# MPI slot used by the master process

n=mpi.universe.size() - 1
print(paste("universe is: ", n))
mpi.spawn.Rslaves(nslaves=n)

# Tell all slaves to return a message identifying themselves
mpi.remote.exec(
    paste(
      "I am",mpi.comm.rank(),
      "of", mpi.comm.size()-1,
      "on host", Sys.info()[c("nodename")]
    )
)
# Tell all slaves to close down, and exit the program
mpi.close.Rslaves()
mpi.quit()
```

will produce output like this:

```
gizmof369[~/tutorial]: mpirun -n 1 rmpi.new.R 
[1] "universe is:  2"
    2 slaves are spawned successfully. 0 failed.
master (rank 0, comm 1) of size 3 is running on: gizmof369 
slave1 (rank 1, comm 1) of size 3 is running on: gizmof370 
slave2 (rank 2, comm 1) of size 3 is running on: gizmof371 
$slave1
[1] "I am 1 of 2 on host gizmof370"

$slave2
[1] "I am 2 of 2 on host gizmof371"

[1] 1
```
