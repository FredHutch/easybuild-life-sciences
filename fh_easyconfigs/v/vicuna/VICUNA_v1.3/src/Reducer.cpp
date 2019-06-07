//============================================================================
// Project     : Diversifier
// Name        : Reducer.cpp
// Author      : Xiao Yang
// Created on  : Aug 17, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#include "Reducer.h"

void Reducer::run (const ReadIndexer& rIndexer, int batchSize,
		const Trimmer& myTrimmer) {

	std::cout << "Reducing ...\n\n";
	std::cout << "\tSuper-shingling...\n\n";
	super_shingling (rIndexer, batchSize, myTrimmer);

	// debugging
	std::cout << "\t\t super-shingled "<< _mySS.size() << " reads ("
			<< (100.0 * _mySS.size())/(2*rIndexer.size()) << " %total)\n\n";

	std::cout << "\tClustering...\n\n";
	clustering (rIndexer.size(), rIndexer);

	//std::cout << "\tValidation...\n\n";
}

// for debugging purposes, introduce rIndexer as the parameter
void Reducer::clustering(int numReads, const ReadIndexer& rIndexer) {

	//1. generate super-shingle-> rIDs structure

	std::cout << "\t\tsuper-shingles --> rIDs\n";

	std::map<uint32_t, ivec_t> ss2rIDs;
	std::map<uint32_t, ivec_t>::iterator it;
	for (int i = 0; i < (int) _mySS.size(); ++ i) {
		if (_mySS[i].rID == -1) continue;
		it = ss2rIDs.find(_mySS[i].fwdSS);
		if (it != ss2rIDs.end())	it->second.push_back(_mySS[i].rID);
		else ss2rIDs[_mySS[i].fwdSS] = ivec_t (1, _mySS[i].rID);

		it = ss2rIDs.find(_mySS[i].rvSS);
		if (it != ss2rIDs.end())	it->second.push_back(_mySS[i].rID);
		else ss2rIDs[_mySS[i].rvSS] = ivec_t (1, _mySS[i].rID);
	}


	//------------debugging -- print out histogram here
	/*
	ivec_t histo;
	int debug_counter = 0;
	for (it = ss2rIDs.begin(); it!= ss2rIDs.end(); ++ it) {
		histo.push_back(it->second.size());
		/*
		if (it->second.size() > 1) {
			std::vector<read_t> reads;
			rIndexer.getReads(std::back_inserter(reads), it->second.begin(), it->second.end());
			std::cout << "cls:" << debug_counter << "\t" << reads.size() << "\t" << it->first << "\n";
			for (int i = 0; i < reads.size(); ++ i) {
				std::cout << reads[i].read << "\n";
				std::cout << xny::get_rvc_str(reads[i].read) << "\n\n";
			}
			debug_counter ++ ;
		}
		*/
	/*
	}
	std::cout << "\n\n";
	std::cout << "\thistogram of ss->rID size\n\n";
	histogram(histo,  (&std::cout));
	std::cout << "\n\t==========\n\n";
	*/ // -----------end debugging -------------------------

	// consolidate the map structure to iivec_t; since ss matters no more
	iivec_t cls_rIDs;
	cls_rIDs.reserve(ss2rIDs.size());
	for (it = ss2rIDs.begin(); it != ss2rIDs.end(); ++ it)
		cls_rIDs.push_back(it->second);
	ss2rIDs.clear();

	// 2. uf structure for clustering.
	// Treat each super-shingle->rIDs as a cluster,
	// where super-shingle is the cluster ID, if a read ID belongs to
	// different clusters, then merge would happen.
	// Use uf structure to cluster super-shingles

	std::cout << "\t\t UF: generate cls of rIDs\n";

	ivec_t ufvec (cls_rIDs.size(), 0); //
	// initialize
	for (int i = 0; i < (int) ufvec.size(); ++ i) ufvec[i] = i;
	// generate rID-> clsIdx (es)
	iivec_t rID2clsIdx (numReads, ivec_t());
	for (int i = 0; i < (int) cls_rIDs.size(); ++ i) {// check every cluster
		for (int j = 0; j < (int) cls_rIDs[i].size(); ++ j)
			rID2clsIdx[cls_rIDs[i][j]].push_back(i);
	}
	// uf
	for (int i = 0; i < (int) rID2clsIdx.size(); ++ i) { // for every rID
		for (int j = 0; j < (int) rID2clsIdx[i].size() - 1; ++ j) {
			int clsIdx1 = rID2clsIdx[i][j];
			int root_1 = find (clsIdx1, ufvec);
			for (int l = j + 1; l < (int) rID2clsIdx[i].size(); ++ l) {
				int clsIdx2 = rID2clsIdx[i][l];
				int root_2 = find (clsIdx2, ufvec);
				if (root_1 != root_2) ufvec[root_2] = root_1;
			}
		}
	}
	// root -> clsIdx(es)
	std::cout << "\t\t Forming clusters\n";

	std::map<int, ivec_t> root2clsIdx;
	std::map<int, ivec_t>::iterator r2c_it;

	for (int i = 0; i < (int) ufvec.size(); ++ i) {
		int root = clsfind (i, ufvec);
		r2c_it = root2clsIdx.find(root);
		if (r2c_it != root2clsIdx.end()) {
			r2c_it->second.push_back(i);
		} else {
			root2clsIdx[root] = ivec_t (1, i);
		}
	}
	// merge rIDs to create final clusters
	std::cout << "\t\tMerging rIDs to form final clusters\n";

	std::vector<iset_t> union_cls_rIDs;
	for (r2c_it = root2clsIdx.begin(); r2c_it != root2clsIdx.end(); ++ r2c_it){
		iset_t union_rIDs;
		for (int i = 0; i < (int) r2c_it->second.size(); ++ i) {
			int clsIdx = r2c_it->second[i];
			union_rIDs.insert(cls_rIDs[clsIdx].begin(), cls_rIDs[clsIdx].end());
		}
		union_cls_rIDs.push_back(union_rIDs);
	}

	//------------debugging -- print out histogram here
	/*
	std::cout << "\t\tdebug: generate histograms\n\n";
	// generate histogram

	ivec_t hist (union_cls_rIDs.size(), 0);
	int totalrIDs = 0;
	for (int i = 0; i < (int) hist.size(); ++ i) {
		hist [i] = union_cls_rIDs[i].size();
		totalrIDs += hist[i];
/*
		// debug print out
		if (union_cls_rIDs[i].size() > 1) {
			std::vector<read_t> reads;
			rIndexer.getReads(std::back_inserter(reads),
					union_cls_rIDs[i].begin(), union_cls_rIDs[i].end());
			std::cout << "cls:" << debug_counter << "\t" << reads.size() << "\n";
			for (int i = 0; i < reads.size(); ++ i) {
				std::cout << reads[i].read << "\n";
				std::cout << xny::get_rvc_str(reads[i].read) << "\n\n";
			}
			debug_counter ++ ;
		}
*/
	/*
	}
	std::cout << "\tSingletons after super-shingling: "
			<< numReads - totalrIDs << "\n";
	std::cout << "\tcls_size\t#clusters\n";
	histogram (hist, (&std::cout));
	*/
}

/* generate super shingles for both fwd and rv strands of each read if
 * the read length is sufficient long
 */

void Reducer::super_shingling (const ReadIndexer& rIndexer, int batchSize,
		const Trimmer& myTrimmer){

	// process reads in batches
	std::vector<read_t> reads;
	int counter = 0;
	ivec_t rIDs;

	for (int i = 0; i < (int) rIndexer.size(); ++ i) {
		rIDs.push_back(i);
		counter ++;
		if (counter >= batchSize) {
			rIndexer.getReads(std::back_inserter(reads), rIDs.begin(),
					rIDs.end());
			// process this batch
			batch_ss (reads, myTrimmer);

			counter = 0;
			reads.clear();
			rIDs.clear();
		}
	}
	if (counter > 0) {
		rIndexer.getReads(std::back_inserter(reads), rIDs.begin(),
				rIDs.end());
		// process the remaining
		batch_ss (reads, myTrimmer);
	}
}


void Reducer::batch_ss (const std::vector<read_t>& reads,
		const Trimmer& myTrimmer) {


	int preSize = reads.size();
	std::vector<ss_t> localSS (preSize, ss_t());
	std::string myRead;
#pragma omp parallel for shared (localSS) private (myRead)
	for (int i = 0; i < preSize; ++ i) {
		myTrimmer.trim (myRead, reads[i]);
		localSS[i].rID = -1;
		if (gen_ss(localSS[i].fwdSS, myRead) &&
				gen_ss(localSS[i].rvSS, xny::get_rvc_str(myRead))) {
			localSS[i].rID = reads[i].ID;
		}
	}
	_mySS.insert(_mySS.end(), localSS.begin(), localSS.end());
}

bool Reducer::gen_ss(uint32_t& ss, const std::string& read) {

	jaz::to_string ts;

	// stage 1
	int l = read.size();

	if (l - _w1 + 1 <= 0) return false;

	std::vector<uint32_t> sh1(l - _w1 + 1, 0);

	// jaz::murmur264 hash
	jaz::djb32 hash;

	for (int i = 0; i < l - _w1 + 1; ++i) {
		sh1[i] = hash(std::string(read.begin() + i, read.begin() + i + _w1));
	}

	std::sort(sh1.begin(), sh1.end());

	// stage 2
	l = sh1.size();
	std::vector<uint32_t> sh2(l - _w2 + 1, 0);

	if (l - _w2 + 1 <= 0) return false;

	/*
	for (int i = 0; i < l - _w2 + 1; ++i) {
		std::string v = ts(sh1[i]);
		for (int j = 1; j < _w2; ++j) {
			v += ' ' + ts(sh1[i + j]);
		}
		sh2[i] = hash(v);
	}*/
    unsigned int v_sz = _w2 * sizeof(uint32_t);

    for (unsigned int i = 0; i < l - _w2 + 1; ++i) {
    		const char* v = reinterpret_cast<const char*>(&sh1[i]);
    		sh2[i] = hash(v, v_sz);
    }


	ss = *std::min_element(sh2.begin(), sh2.end());

	return true;
}
