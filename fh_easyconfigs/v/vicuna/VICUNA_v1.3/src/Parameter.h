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
#include "xutil.h"
#include "KmerType.h"


class Parameter{

public:
	//-------------- Trimmer ----------------------
	std::string t_iVecFileName, t_ilogName;
	int t_minMSz;
	int t_minIMSz;
	int t_maxOverhang;
	int t_minReadSz;

	//-------------- Profiler ----------------------
	std::string p_iMSAFileName;
	int p_binNum;
	int p_K, p_HD, p_blockNum, p_minSpan;
	std::string p_oRMapFileName;

	//-------------- Contiger ----------------------

	// word size for the first and the second shingling
	int c_w1, c_w2;
	// maxmum % of divergence between read & consensus
	int c_percent_diverge;
	// number of read end bps to be ignored, due to insufficient trimming,
	// PCR artifacts, sequencing errors etc.
	int c_read_overhang;

	// to determine the reliable region of a consensus
	int c_min_profile_col_weight;
	int c_min_consensus_base_ratio;
	int c_max_contig_overhang;
	// merging criteria of two contigs
	int c_seed_kmer_len;
	int c_min_contig_overlap;
	int	c_min_contig_links;
	int c_min_identity;
	int c_min_perc_pol;
	int c_max_variant_len;
	//-------------- General Assembly --------------
	std::string pFqDir, npFqDir, pFaDir, npFaDir, oDIRNm;
	int min_output_contig_len;
	int libSzLb, libSzUb;
	int 	batchSize;

	Parameter(int argc, char** argv): argnum(argc), arg(argv){
		if (argnum != 2) printUsage(argv[0]);
         std::ifstream input(arg[1]);
        if (!input) { std::cout << "cannot open " << arg[1] << "\n"; exit(1);}
        initialize();
        getPara (input);
		printSpec();
		input.close();
	}

private:
	int argnum;
	char** arg;
	void printUsage(char* exe) {
			std::cout << "\n--------------------------------------------------------\n";
			std::cout << "Usage: " << exe << " [Config File]\n";
			std::cout << "----------------------------------------------------------\n\n";
			exit(1);
	}

	void initialize () {
		// set default values
        // Trimmer
        t_iVecFileName = ""; t_ilogName = "";
        t_minMSz = 7; t_minIMSz = 15;
    		t_maxOverhang = 4; t_minReadSz = 25;
    		// Contiger
    		c_w1 = 12;
    		c_w2 = 5;
    		c_percent_diverge = 10;
    		c_read_overhang = 4;
    		c_min_profile_col_weight = 5;
    		c_min_consensus_base_ratio = 85;
    		c_max_contig_overhang = 10;
    		c_seed_kmer_len = 12;
    		c_min_contig_overlap = 25;
    		c_min_contig_links = 3;
    		c_min_identity = 90;
    		c_min_perc_pol = 5;
    		c_max_variant_len = 20;
		// Profiler
    		p_binNum = 20; p_K = 15;
    		p_blockNum = 5; p_HD = 1;
    		p_minSpan = 75;
    		p_iMSAFileName = "";
    		p_oRMapFileName = "";
		// Assembly
		libSzLb = libSzUb = -1;
		oDIRNm = "./";
		pFqDir = npFqDir = pFaDir = npFaDir = "";
		min_output_contig_len = 300;
		batchSize = 2000000;
	}

	void getPara (std::ifstream& input) {
	   std::string line, s1;
	   std::istringstream buf;
		while (std::getline(input, line)) {
			buf.clear();
			buf.str(line);

			if (buf >> s1) {
				if (s1 == "vectorFileName") { //Trimmer
					buf >> t_iVecFileName;
				} else if (s1 == "trimLogFileName"){
					buf >> t_ilogName;
				} else if (s1 == "minMSize") {
					buf >> t_minMSz;
				} else if (s1 == "minInternalMSize") {
					buf >> t_minIMSz;
				} else if (s1 == "maxOverhangSize") {
					buf >> t_maxOverhang;
				} else if (s1 == "minReadSize") {
					buf >> 	t_minReadSz;
				}
				else if (s1 == "w1") {		// Contiger
					buf >> c_w1;
				} else if (s1 == "w2") {
					buf >> c_w2;
				} else if (s1 == "Divergence") {
		    			buf >> c_percent_diverge;
				} else if (s1 == "min_profile_col_weight") {
		    			buf >> c_min_profile_col_weight;
				} else if (s1 == "min_consensus_base_ratio") {
		    			buf >> c_min_consensus_base_ratio;
				} else if (s1 == "max_contig_overhang") {
		    			buf >> c_max_contig_overhang;
		    		} else if (s1 == "max_read_overhang") {
		    			buf >> c_read_overhang;
				} else if (s1 == "seed_kmer_len") {
					buf >> c_seed_kmer_len;
					if (c_seed_kmer_len > 16 || c_seed_kmer_len < 9) {
						abording("set seed_kmer_len range in [9, 16]");
					}
				} else if (s1 == "min_contig_overlap") {
					buf >> c_min_contig_overlap;
				} else if (s1 == "min_contig_links") {
				    buf >> c_min_contig_links;
				} else if (s1 == "min_identify") {
					buf >> c_min_identity;
				} else if (s1 == "min_perc_polymorphism"){
					buf >> c_min_perc_pol;
				} else if (s1 == "max_variant_len"){
					buf >> c_max_variant_len;
				}else if (s1 == "MSAFileName") { //Profiler
					buf >> p_iMSAFileName;
				} else if (s1 == "binNumber") {
					buf >> p_binNum;
				} else if (s1 == "kmerLength") {
					buf >> p_K;
				} else if (s1 == "maxHD") {
					buf >> p_HD;
				} else if (s1 == "p_minSpan") {
					buf >> p_minSpan;
				} else if (s1 == "blockNumber") {
					buf >> p_blockNum;
				} else if (s1 == "rMapFileName") {
					buf >> p_oRMapFileName;
				}
				else if (s1 == "pFqDir") { //general assembly
					buf >> pFqDir;
				} else if (s1 == "npFqDir" ) {
					buf >> npFqDir;
				} else if (s1 == "pFaDir") {
					buf >> pFaDir;
				} else if (s1 == "npFaDir") {
					buf >> npFaDir;
				} else if (s1 == "batchSize") {
					buf >> batchSize;
					if (batchSize % 2 != 0) ++ batchSize;
				} else if (s1 == "LibSizeLowerBound") {
					buf >> libSzLb;
				} else if (s1 == "LibSizeUpperBound") {
					buf >> libSzUb;
				} else if (s1 == "min_output_contig_len") {
					buf >> min_output_contig_len;
				} else if (s1 == "outputDIR") {
					buf >> oDIRNm;
				}
			} // if
		} // while
	} // getPara

	void printSpec() {
		std::cout << "\n--------------------------------------------------------\n";
		std::cout << "Program runs with the following Parameter setting:\n\n";
		std::cout << "\t===== Trimmer =====\n\n";
		std::cout << "\tvectorFileName\t" << t_iVecFileName << "\n";
		std::cout << "\ttrimLogFileName\t" << t_ilogName << "\n";
		std::cout << "\tminMSize\t" << t_minMSz << "\n";
		std::cout << "\tminInternalMSize\t" << t_minIMSz << "\n";
		std::cout << "\tmaxOverhangSize\t" << t_maxOverhang << "\n";
		std::cout << "\tminReadSize\t" << t_minReadSz << "\n";

		std::cout << "\n\t===== Profiler =====\n\n";
		std::cout << "\tMSAFileName\t" << p_iMSAFileName << "\n";
		std::cout << "\tbinNumber\t" << p_binNum << "\n";
		std::cout << "\tkmerLength\t" << p_K << " (encode using "
				<< sizeof(kmer_t) << " bytes)\n";
		std::cout << "\tmaxHD\t" << p_HD << "\n";
		std::cout << "\tminSpan\t" << p_minSpan << "\n";
		std::cout << "\tblockNumber\t" << p_blockNum << "\n";
		std::cout << "\trMapFileName\t" << p_oRMapFileName << "\n";

		std::cout << "\n\t===== Contiger =====\n\n";
		std::cout << "\tw1\t" << c_w1 << "\n";
		std::cout << "\tw2\t" << c_w2 << "\n";
		std::cout << "\tDivergence\t" << c_percent_diverge << "\n";
		std::cout << "\tmax_read_overhang\t" << c_read_overhang << "\n";
		std::cout << "\tmin_profile_col_weight\t" << c_min_profile_col_weight << "\n";
		std::cout << "\tmin_consensus_base_ratio\t" << c_min_consensus_base_ratio << "\n";
		std::cout << "\tmax_contig_overhang\t" << c_max_contig_overhang << "\n";
		std::cout << "\tseed_kmer_len\t" << c_seed_kmer_len << "\n";
		std::cout << "\tmin_contig_overlap\t" << c_min_contig_overlap << "\n";
		std::cout << "\tmin_contig_links\t" << c_min_contig_links << "\n";
		std::cout << "\tmin_identity\t" << c_min_identity << "\n";
		std::cout << "\tmin_perc_polymorphism\t" << c_min_perc_pol << "\n";
		std::cout << "\tmax_variant_len\t" << c_max_variant_len << "\n";

		std::cout << "\n\t===== Assembly =====\n\n";

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
		std::cout << "\tbatchSize\t" << batchSize << "\n";
		std::cout << "\tLibSizeLowerBound\t" << libSzLb << "\n";
		std::cout << "\tLibSizeUpperBound\t" << libSzUb << "\n";
		std::cout << "\tmin_output_contig_len\t" << min_output_contig_len << "\n";
		std::cout << "\toutputDIR\t" << oDIRNm << "\n";
		std::cout << "\n--------------------------------------------------------\n\n";
	}

};


#endif /* PARAMETER_H_ */
