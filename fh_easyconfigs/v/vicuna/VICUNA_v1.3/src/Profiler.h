//============================================================================
// Project     : Diversifier
// Name        : Profiler.h
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef PROFILER_H_
#define PROFILER_H_

#include <iterator>
#include "Parameter.h"
#include "xutil.h"
#include "xny/seq_manip.hpp"
#include "KmerType.h"
#include "Trimmer.h"
#include "ReadIndexer.h"
/*
 *
 */

typedef struct BINSHIP {
	BINSHIP (int id, const std::vector<u8pair_t>& vec):
		rID (id), vec_bin_weight (vec) {};
	BINSHIP () {}
	int rID;
	std::vector<u8pair_t> vec_bin_weight;
} bin_t;

class Profiler {
public:
	Profiler(const Parameter& myPara) {
		iMSAFileName_ = myPara.p_iMSAFileName;

		binNum_ = myPara.p_binNum;
		K_ = myPara.p_K;
		HD_ = myPara.p_HD;
		blockNum_ = myPara.p_blockNum;
		minSpan_ = myPara.p_minSpan;
		oRMapFileName_ = myPara.p_oRMapFileName;
		batchSize_ = myPara.batchSize;
		libSzUb_ = myPara.libSzUb;
	}
	virtual ~Profiler();
	void run (const ReadIndexer& rIndexer, const Trimmer& myTrimmer);

	void get_binned_rIDs (ivec_t& rIDs);

private:

	std::vector<ipair_t> _initCoord; // initial boundaries for dividing bins

	std::vector<std::pair<kmer_t, uint8_t> > _kmers; // stores sorted (kmerID, binID) for each bin
	std::vector<kmer_t> _masks; // masks used for sorting k-spectrum
						// kmer indices are divided into c blocks,
	                    // d blocks a masked, and k-spectrum are sorted
						// w.r.t. the remaining indices
	std::vector<std::pair<kmer_t, uint8_t> > _duplTb; // stores all sorted
	                    // k-spectrums

	//std::vector<std::vector<ipair_t> > _binship; // _binship[i]:
					  //stores all  (bin, weight) pairs read_i may belong to;
					  // After resolve ambiguity, a unique pair is assigned
	std::vector<bin_t> binship_;

	/* input parameters */
	std::string iMSAFileName_, oRMapFileName_;
	int binNum_, K_, HD_, blockNum_, minSpan_, batchSize_, libSzUb_;

	void readMSA(strvec_t& seqs);

	void profiling (const std::vector<std::string>& seqs);
	void initCoord (const strvec_t& seqs);
	void identifyBoundary (std::vector<std::vector<ipair_t> >& coord,
			const std::vector<std::string>& seqs);
	std::vector<ipair_t> adjust (std::vector<ipair_t> coord,
								const std::string& seq);
	void profileDupl (); // duplicate profile to facilitate binary search

	void rslvAmbig (int libSizeub, const ReadIndexer& rIndexer);
	void readAssignment (const ReadIndexer& rIndexer,
						const Trimmer& myTrimmer);
	void mapping (const std::vector<read_t>& reads,
				  const Trimmer& myTrimmer);
	void outputBinnedRIDs(const std::string& oFileNm);

};


template <typename int_t>
struct compKT{
	int_t mask;
	compKT (int_t m): mask (m) {};

	bool operator() (const std::pair<int_t, uint8_t>& lh,
			const std::pair<int_t, uint8_t>& rh) const {
		return ((lh.first & mask) < (rh.first & mask));
	}
};

template <typename int_t>
struct compKP{
	bool operator() (const std::pair<int_t, uint8_t>& lh,
			const std::pair<int_t, uint8_t>& rh) const {
		return (lh.first < rh.first);
	}
};
#endif /* PROFILER_H_ */
