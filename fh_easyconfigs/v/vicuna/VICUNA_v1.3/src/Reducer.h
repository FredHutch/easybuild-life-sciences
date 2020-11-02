//============================================================================
// Project     : Diversifier
// Name        : Reducer.h
// Author      : Xiao Yang
// Created on  : Aug 17, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef REDUCER_H_
#define REDUCER_H_

/*
 *
 */

#include "ReadIndexer.h"
#include "Trimmer.h"
#include "jaz/functional_add.hpp"
#include "jaz/string_add.hpp"

typedef struct ss{
	ss (int id, uint32_t fwd, uint32_t rv):
		rID(id), fwdSS(fwd), rvSS(rv){}
	ss (){};
	int rID;
	uint32_t fwdSS;
	uint32_t rvSS;
} ss_t;

//inline bool operator<(const ss_t& lhs, const ss_t& rhs) {
//    return lhs.rID < rhs.rID;
//} // operator<


class Reducer {
public:
	Reducer(int w1, int w2): _w1(w1), _w2 (w2) {}
	virtual ~Reducer(){};
	void run (const ReadIndexer& rIndexer, int batchSize,
			const Trimmer& myTrimmer);
private:
	int _w1, _w2;
	// read without sufficient length will be ignored
	std::vector<ss_t> _mySS;

	void super_shingling (const ReadIndexer& rIndexer, int batchSize,
			const Trimmer& myTrimmer);
	void clustering (int numReads, const ReadIndexer& rIndexer);

	void batch_ss (const std::vector<read_t>& reads, const Trimmer& myTrimmer);

	bool gen_ss(uint32_t& ss, const std::string& read);
};

#endif /* REDUCER_H_ */
