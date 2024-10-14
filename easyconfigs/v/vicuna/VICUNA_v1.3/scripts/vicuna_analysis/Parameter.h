//============================================================================
// Project     : Diversifier
// Name        : Parameter.h
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>

/* Contig regions of interest */
struct region_t {
	region_t (int id, int s, int e): cID (id), start (s), end (e) {}
	int cID;
	int start;
	int end;
};

class Parameter{

public:
	/* [inDirNm]: contains input read files
	 * [trimlogNm]: trim logs
	 * [alnNm]: alignment file
	 * [outDirNm]: output dir
	 */
	std::string trimlogNm, alnNm, oDir;
	std::string pFqDir, npFqDir, pFaDir, npFaDir;
	int num_region;
	std::vector<region_t> regions;

	int lfv_max_len, lfv_freq;

	Parameter(int argc, char** argv): argnum(argc), arg(argv){
		if (argnum != 2) printUsage(argv[0]);
         std::ifstream input(arg[1]);
        if (!input) {
        		std::cout << "cannot open " << arg[1] << "\n";
        		exit(1);
        }
        initialize();
        getPara (input);
		printSpec();

		input.close();
	}

private:
	int argnum;
	char** arg;
	void printUsage(char* exe) {
		std::cout << "\n---------------------------------------------\n";
		std::cout << "Usage: " << exe << " [Config File]\n";
		std::cout << "-----------------------------------------------\n\n";
		exit(1);
	} // printUsage

	void initialize () {
		trimlogNm = alnNm = oDir = "";
		pFqDir = npFqDir = pFaDir = npFaDir = "";
		num_region = 0;
		lfv_freq = 5;
		lfv_max_len = 20;
	} // initialize

	void getPara (std::ifstream& input) {
	   std::string line, s1;
	   std::istringstream buf;
		while (std::getline(input, line)) {
			buf.clear();
			buf.str(line);

			if (buf >> s1) {
				if (s1 == "trim_log_file") { //Trimmer
					buf >> trimlogNm;
				} else if (s1 == "aln_file"){
					buf >> alnNm;
				} else if (s1 == "outputDIR") {
					buf >> oDir;
				} else if (s1 == "pFqDir") { //general assembly
					buf >> pFqDir;
				} else if (s1 == "npFqDir" ) {
					buf >> npFqDir;
				} else if (s1 == "pFaDir") {
					buf >> pFaDir;
				} else if (s1 == "npFaDir") {
					buf >> npFaDir;
				} else if (s1 == "num_region") {
					buf >> num_region;
					for (int i = 0; i < num_region; ++ i) {
						int cID, start, end;
						input >> cID >> start >> end;
						regions.push_back(region_t (cID, start, end));
					}
				} else if (s1 == "lfv_freq") {
					buf >> lfv_freq;
				} else if (s1 == "lfv_max_len") {
					buf >> lfv_max_len;
				}
			}
		} // while
		if (trimlogNm.empty()) {
			std::cout << "trim_log_file is not specified, exit.\n";
			exit(1);
		}
		if (alnNm.empty()) {
			std::cout << "aln_file is not specified, exit.\n";
			exit(1);
		}
		if (oDir.empty()) {
			std::cout << "outputDIR is not specified, exit.\n";
			exit(1);
		}
		if (pFqDir.empty() && npFqDir.empty()
				&& pFaDir.empty() && npFaDir.empty()) {
			std::cout << "no input read file specified, exit.\n";
			exit(1);
		}
	} // getPara

	void printSpec() {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Program runs with the following Parameter setting:\n\n";
		if (!pFqDir.empty()) {
			std::cout << "\tpFqDir\t" << pFqDir << "\n";
		}
		if (!npFqDir.empty()) {
			std::cout << "\tnpFqDir\t" << npFqDir << "\n";
		}
		if (!pFaDir.empty()) {
			std::cout << "\tpFaDir\t" << pFaDir << "\n";
		}
		if (!npFaDir.empty()) {
			std::cout << "\tnpFaDir\t" << npFaDir << "\n";
		}
		std::cout << "\ttrim_log_file\t" << trimlogNm << "\n";
		std::cout << "\taln_file\t" << alnNm << "\n";
		std::cout << "\toutputDIR\t" << oDir << "\n";
		if (num_region > 0) {
			std::cout << "\tnum_region\t" << num_region << "\n";
			std::cout << "\t\tcID\tStart\tEnd\n";
			for (int i = 0; i < num_region; ++ i) {
				std::cout << "\t\t" << regions[i].cID << "\t"
						<< regions[i].start << "\t"
						<< regions[i].end << "\n";
			}
		}
		std::cout << "\tlfv_freq\t" << lfv_freq << "\n";
		std::cout << "\tlfv_max_len\t" << lfv_max_len << "\n";

		std::cout << "\n--------------------------------------------------------\n\n";
	} // printSpec

};


#endif /* PARAMETER_H_ */
