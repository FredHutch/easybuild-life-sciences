//============================================================================
// Project     : DivAnalysis
// Name        : contig.h
// Author      : Xiao Yang
// Created on  : Dec 13, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef CONTIG_MANIP_H_
#define CONTIG_MANIP_H_

#include <iomanip>
#include <vector>
#include <set>
#include "jaz/string_add.hpp"
#include "xny/seq_manip.hpp"
#include "readEntry.h"
#include "ReadIndexer.h"
#include "Parameter.h"

typedef std::vector<int> ivec_t;
typedef std::set<int> iset_t;
typedef std::pair<int, int> ipair_t;

typedef struct ContigRead {
	ContigRead (){}
	std::string rName;
	int rID;			    // read ID
	int rLen;
	int dist_to_beg; 	// read start position w.r.t. contig start
	bool is_fwd_strand; // w.r.t. the input strand
	bool is_paired;
	std::vector<ipair_t> indels; // indel: <start, len>
	int indel_length () const {
		int indel_len = 0;
		for (int i = 0; i < (int) indels.size(); ++ i) {
			indel_len += indels[i].second;
		}
		return indel_len;
	}
} cread_t;


typedef struct ContigAln {
	ContigAln(){}
	int cID;
	int cLen;
	std::vector<cread_t> rElems;
} contig_aln_t;

class ContigManip {
public:
	ContigManip (const Parameter& myPara){
		_oDir = myPara.oDir;
		_aln_fNm = myPara.alnNm;
		_num_region = myPara.num_region;
		if (_num_region > 0) _regions = myPara.regions;
		_lfv_max_len = myPara.lfv_max_len;
		_lfv_freq = myPara.lfv_freq;

		readinput();
	}
	void output_contigs (const ReadIndexer& rIndexer,
			const Trimmer& myTrimmer);

	void output (const ReadIndexer& rIndexer, const Trimmer& myTrimmer);

private:
	std::vector<contig_aln_t> _cdb;

	/* input parameters */
	std::string _oDir;
	std::string _aln_fNm;
	int _num_region;
	std::vector<region_t> _regions;
	int _lfv_max_len, _lfv_freq;

	void readinput ();
	void id_lfv (std::vector<ipair_t>& lfv, const ivec_t& profile,
			int cID, int cLen);
	void get_lfv_free_contig (std::string& rcons,
			const std::string& consensus, const std::vector<ipair_t>& lfv);

	/* print out alignment; if [start] && [end] != -1, then print out
	 * read layout between position [start] and [end] wrt contig [cID] */
	void print_aln(std::ofstream& oHandle,
			int cID, const ReadIndexer& rIndexer,
			const Trimmer& myTrimmer, bool to_print_profile,
			bool to_print_reads, int start = -1, int end = -1);

	/* profile */
	void generate_profile (std::string& consensus,
		ivec_t& profile, int cID, const std::vector<rentry_t>& reads);
	void update_profile (std::string& consensus, ivec_t& profile,
	  cread_t cread, const std::string& read_seq, int add);
	void update_profile_column (std::string& consensus,
			ivec_t& profile, int col, char newchar, int value);
	void print_profile(std::ofstream& oHandle,
			const std::string& consensus, const ivec_t& profile);

	/* reads */
	void print_reads (std::ofstream& oHandle, int cID,
			const std::vector<rentry_t>& reads, int start, int end);

	/* get all rIDs from [cID] */
	void get_rIDs (iset_t& rIDs, int cID) {
		if (cID < 0 || cID >= (int) _cdb.size()) {
			std::cout << "err: get_rIDs -- cID = " << cID
					<< " is out of range, exiting...\n";
			exit(1);
		}
		int sz = _cdb[cID].rElems.size();

		for (int i = 0; i < sz; ++ i) {
			rIDs.insert(_cdb[cID].rElems[i].rID);
		}
	}

	void get_rIDs (ivec_t& rIDs, int cID) {
		if (cID < 0 || cID >= (int) _cdb.size()) {
			std::cout << "err: get_rIDs -- cID = " << cID
					<< " is out of range, exiting...\n";
			exit(1);
		}
		int sz = _cdb[cID].rElems.size();
		rIDs.resize(sz);
		for (int i = 0; i < sz; ++ i) {
			rIDs[i] = _cdb[cID].rElems[i].rID;
		}
	}
};

#endif /* CONTIG_MANIP_H_ */
