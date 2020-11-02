//============================================================================
// Project     : DivAnalysis
// Name        : main.cpp
// Author      : Xiao Yang
// Created on  : Nov 12, 2011
// Version     :
// Copyright   :
// Description :	Given paired-fasta/fastq files as input, randomly sample
//				x number of read-pairs and generate two new paired files
//============================================================================


#include <iostream>
#include "Parameter.h"
#include "xny/file_manip.hpp"
#include "ReadIndexer.h"
#include "contig_manip.h"

int main (int argc, char** argv){

	Parameter myPara(argc, argv);

	/* getting trimming information */
	Trimmer myTrimmer (myPara.trimlogNm);

	/* indexing reads */
	ReadIndexer rIndexer (myPara.pFqDir, myPara.npFqDir,
			myPara.pFaDir, myPara.npFaDir);

	/* reading contig alignment information */
	ContigManip cManiper (myPara);
	cManiper.output(rIndexer, myTrimmer);
	cManiper.output_contigs (rIndexer, myTrimmer);
	std::cout << "DONE!\n\n";
	return (EXIT_SUCCESS);
}


