## PEER

(PEER) [https://www.sanger.ac.uk/science/tools/peer]

PEER is a collection of libraries that comes from a single repository. Within a subdirectory
of the repository is an R library.

To build both the library, and the R library mulitble steps are required.

Build peer-mastar first.  This will create the libraries and install the R library files.
Create a tarball from the R directory of thd of the install location of peer-master.
Use the 'peer' for the tarball and R libray

 - tar -czC %(installdir)/peer-master/R -f peer-1.0.0.tar.gz peer

Build the peer-1.0.0  RPackage using peer-master as a dependency.
