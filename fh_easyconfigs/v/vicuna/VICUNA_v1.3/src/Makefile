MYPATH=/Users/xiaoyang/Desktop/XIAO/program/LIBRARY/ncbi_cxx--7_0_0/ncbi_install/
COMPILER=/usr/bin/g++

SEQ_LIBS = -lseq -lseqcode -lsequtil
SOBJMGR_LDEP = -lgenome_collection -lseqset -lseqedit $(SEQ_LIBS) -lpub -lmedline -lbiblio -lgeneral -lxser -lxutil -lxncbi
SOBJMGR_LIBS = -lxobjmgr $(SOBJMGR_LDEP)
LIB = -lxalgoalignnw -ltables $(SOBJMGR_LIBS)

all:
	$(COMPILER) -fopenmp -O3 -I$(MYPATH)/include/ncbi-tools++ -L$(MYPATH)/lib *.cpp -o ../bin/vicuna-omp-v1.0 -lpthread $(LIB)
