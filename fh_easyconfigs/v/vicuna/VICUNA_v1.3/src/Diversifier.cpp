//============================================================================
// Project     : Diversifier
// Name        : Diversifier.cpp
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include "Parameter.h"
#include "Trimmer.h"
#include "Profiler.h"
#include "ReadIndexer.h"
#include "xny/file_manip.hpp"
#include "Contiger.h"

int main (int argc, char** argv){

	double timing = get_time();
	double start_time = timing;
	Parameter myPara(argc, argv);

	// Indexing
	ReadIndexer rIndexer (myPara.pFqDir, myPara.npFqDir,
						  myPara.pFaDir, myPara.npFaDir);

    print_time("Indexing done !!!\t", timing);
	std::cout << "-----------------------------------------------------------\n\n";

	// Trimming

	Trimmer myTrimmer (myPara);
	myTrimmer.run (rIndexer);

    print_time("Trimming done !!!\t", timing);
	std::cout << "-----------------------------------------------------------\n\n";

	// Profiling
	Profiler myProfiler (myPara);
	myProfiler.run (rIndexer, myTrimmer);

    print_time("Profiling done !!!\t", timing);
	std::cout << "-----------------------------------------------------------\n\n";

	ivec_t rIDs;
	Contiger myContiger (myPara);

	bool is_filter_on = false;
	if (!myPara.p_iMSAFileName.empty()) {
		/* getting mapped read IDs from myProfiler */
		myProfiler.get_binned_rIDs(rIDs);
		is_filter_on = true;
	} else {
		/* no profiling stage, all input reads are considered */
		int total_reads = rIndexer.size();
		rIDs.resize (total_reads);
		for (int i = 0; i < total_reads; ++ i) rIDs[i] = i;
	}

	myContiger.run (is_filter_on, rIDs, myTrimmer, rIndexer);

    print_time("Contiger done!!!\t", timing);
	std::cout << "-----------------------------------------------------------\n\n";

    print_time("Whole program takes \t", start_time);
	std::cout << "DONE!\n";
	return (EXIT_SUCCESS);
}


