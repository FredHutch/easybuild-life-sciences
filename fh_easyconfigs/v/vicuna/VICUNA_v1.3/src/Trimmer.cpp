//============================================================================
// Project     : Diversifier
// Name        : Trimmer.cpp
// Author      : Xiao Yang
// Created on  : Aug 2, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#include "Trimmer.h"
#include "ReadIndexer.h"

Trimmer::~Trimmer() {
}


Trimmer::Trimmer(const Parameter& myPara){
	maxOverhang_ = myPara.t_maxOverhang;
	minInternalMSz_ = myPara.t_minIMSz;
	minMSz_ = myPara.t_minMSz;
	minReadSz_ = myPara.t_minReadSz;
	_logName = myPara.t_ilogName;
	_vecName = myPara.t_iVecFileName;
	_batchsz = myPara.batchSize;
	_oDir = myPara.oDIRNm;

	numHeads_ = 0;
	numTails_ = 0;
	numFull_ = 0;
	numlowComplexity_ = 0;
}

/* @brief	Identify readID in trim log, which is used to trim the read
 * 			when applicable, trimming result is stored in [myRead]
 */
void Trimmer::trim (std::string& myRead, const read_t& read) const {

	myRead = read.read;

	typedef 	std::vector<trimlog_t>::const_iterator iter_t;

	std::pair<iter_t, iter_t> bounds = std::equal_range(
			_logs.begin(), _logs.end(), trimlog_t(read.ID));

	if (bounds.first != bounds.second) { // found
		if (bounds.first->from == 0) {
			if (bounds.first->to == (int) myRead.length() - 1) { // rm read
				myRead = "";
			} else { // trim head
				myRead = myRead.substr(bounds.first->to + 1,
					(int) myRead.length() - (bounds.first->to + 1));
			}
		} else { // trim tail
			myRead = myRead.substr(0, bounds.first->from);
		}
	}
}

void Trimmer::run (const ReadIndexer& rIndexer) {

	std::cout << "Trimming ...\n\n";

	if (!_logName.empty()) {
		read_trim_log ();
		return;
	}

	/* vector trimming */
	strvec_t vectors;

	if (_vecName.empty()) std::cout << "\tno vector trimming applied.\n\n";
	else {
		/* read in vectors */
		readFasta (std::back_inserter(vectors), _vecName);

		/* print out vectors */
		std::cout <<"\tVectors used for trimming: \n";
		for (int i = 0; i < (int) vectors.size(); ++ i) {
			std::cout << "\t" << vectors[i] << "\n";
		}
		std::cout << "\n";
	}

	// 2. process reads in batches
	std::vector<read_t> reads;
	int counter = 0;
	ivec_t rIDs;
	for (int i = 0; i < (int) rIndexer.size(); ++ i) {
		rIDs.push_back(i);
		counter ++;
		if (counter >= _batchsz) {
			rIndexer.getReads(std::back_inserter(reads), rIDs.begin(),
					rIDs.end());

			// process this batch
			vBatchTrim (reads, vectors);
			counter = 0;
			reads.clear();
			rIDs.clear();
		}
	}
	if (counter > 0) {
		rIndexer.getReads(std::back_inserter(reads), rIDs.begin(),
				rIDs.end());
		// process the remaining
		vBatchTrim (reads, vectors);
	}

	// 3. clean out logs (by removing entries with rID = -1)

	std::vector<trimlog_t> tmp_logs;
	for (int i = 0; i < (int) _logs.size(); ++ i) {
		if (_logs[i].readID != -1) tmp_logs.push_back(_logs[i]);
	}
	_logs = tmp_logs;

	// 4. output
	output();

} // run

/***********************************************************************
 * @brief	Read trim logs from an input file in tab delimited format:
 * 			[rID] [start_pos] [end_pos]
 ***********************************************************************/
void Trimmer::read_trim_log () {
	std::ifstream iHandle (_logName.c_str());
	if (!iHandle.good()) {
		std::cout << "\tCan't open trim log file: " + _logName << "\n";
		exit(1);
	} else {
		std::cout << "\tReading trim logs from file: " << _logName << "\n";
		int rID, start, end;
		while (!iHandle.eof()) {
			iHandle >> rID >> start >> end;
			_logs.push_back(trimlog_t(rID, start, end, 0));
		}
	}
} // read_trim_log

/* @brief	1) For vector trimming, trim the ends or discard the whole read:
 * 			look for exact match with seeds, then extend in both DIRs
 * 			2) For low complexity read trimming
 */
void Trimmer::vBatchTrim (const std::vector<read_t>& reads,
		const strvec_t& vectors){

	/* Preparation for vector trimming: for each kmer of each vector,
	 * use [seedInfo] to store [kmerID -> (whichVec, pos, dir)] */
	typedef std::map<uint32_t, std::vector<seed_t> > record_t;
	record_t seedInfo;
	for (int i = 0; i < (int) vectors.size(); ++ i) {
		// generate (minMSz_)-spectrum for fwd & rvc strands of vector[i]
		uvec_t seeds;
		xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t>
			(vectors[i], std::back_inserter(seeds), minMSz_, 2);

		for (int j = 0; j < (int) seeds.size(); ++ j) {
			record_t::iterator it  = seedInfo.find(seeds[j]);
			if (it != seedInfo.end()) {
				if (j % 2 == 0)
					it->second.push_back(seed_t(i, j/2, true));
				else it->second.push_back(seed_t(i, j/2, false));
			} else {
				if (j % 2 == 0) {
					seedInfo[seeds[j]] =
							std::vector<seed_t> (1, seed_t(i, j/2, true));
				}
				else {
					seedInfo[seeds[j]] =
						std::vector<seed_t> (1, seed_t(i, j/2, false));
				}
			}
		} // j
	} // i

	/* parallel trimming and store the result in [local_logs] */

	int size = reads.size();

	std::vector<trimlog_t> local_logs (size, trimlog_t());

	#pragma omp parallel for shared (local_logs)
	for (int i = 0; i < size; ++ i) {

		/* low complexity read trimming */
		if (islowComplexity(reads[i].read)) {

			++ numlowComplexity_;

			local_logs[i] = trimlog_t(reads[i].ID, 0,
					reads[i].read.length() - 1, reads[i].paired);

			continue;
		}

		/* vector trimming */

		if (vectors.size() == 0) continue;

		uvec_t rkmers;
		xny::get_bitkmer<std::back_insert_iterator<uvec_t>, uint32_t>
			(reads[i].read, std::back_inserter(rkmers), minMSz_, 3);

		for (int j = 0; j < (int) rkmers.size(); ++ j) {

			record_t::iterator it  = seedInfo.find(rkmers[j]);

			// check if we can do trimming using current kmer: kmers[j]
			if (it != seedInfo.end()) {

				trimlog_t myLog;
				if (trimmable (myLog, it->second, vectors, reads[i], j)) {
					local_logs[i] = myLog;
					break;
				}
			}

		} // for j

	} // for i

	_logs.insert(_logs.end(), local_logs.begin(), local_logs.end());
} // vBatchTrim


 /*
 * For every kmer shared between read and any vector,
 * return once the first trimming event happens
 */
bool Trimmer::trimmable (trimlog_t& mylog, const std::vector<seed_t>& kinfo,
		const strvec_t& vectors, const read_t& myRead, int rPos) {

	int rBegin = rPos, rEnd = rPos + minMSz_ - 1;
	std::string rSeq = myRead.read;
	int rLen = rSeq.length();

	std::vector<seed_t>::const_iterator it = kinfo.begin();

	for (; it != kinfo.end(); ++ it) { // check for each matching kmer

		std::string myVec = vectors[it->vecID];// vector string
		int vBegin = it->start, vEnd = it->start + minMSz_ - 1;
		bool vDir = it->dir;					  // dir

		std::string rPrf = rSeq.substr(0, rBegin),
				    rSuf = rSeq.substr(rEnd+1, rLen - rEnd),
					vPrf = myVec.substr(0, vBegin),
					vSuf = myVec.substr(vEnd+1, myVec.length() - vEnd);

		if (vDir) { // fwd
			rEnd += xny::common_pref_len (rSuf, vSuf);
			rBegin -= xny::common_suf_len (rPrf, vPrf);
		} else { // rv
			rEnd += xny::common_pref_len (rSuf, xny::get_rvc_str (vPrf));
			rBegin -= xny::common_suf_len (rPrf, xny::get_rvc_str(vSuf));
		}

		// check for trimming criteria
		if (rBegin <= maxOverhang_) {
			if (rLen - rEnd - 1 <= maxOverhang_) { // rm read
				++ numFull_;
				//myRead.read = std::string(1, 'x');
				//log = log_t(myRead.header, rBegin, rEnd, it->dir);
				mylog = trimlog_t(myRead.ID, 0, rLen-1, myRead.paired);
				return true;
			} else {
				if (rLen - rEnd - 1 >= minReadSz_) { // trim prefix
					//rSeq = rSeq.substr(rEnd + 1, std::max(0, rLen - rEnd - 1));
					mylog = trimlog_t (myRead.ID, 0, rEnd, myRead.paired);
					++ numHeads_;
				} else { // rm read
					mylog = trimlog_t (myRead.ID, 0, rLen - 1, myRead.paired);
					++ numFull_;
					//myRead.read = std::string(1, 'x');
				}
				// log = log_t(myRead.header, rBegin, rEnd, it->dir);
				return true;
			} // if ... else
		} else {
			if (rLen - rEnd - 1 <= maxOverhang_) {
				if (rBegin >= minReadSz_) { // trim suffix
					mylog = trimlog_t (myRead.ID, rBegin, rLen - 1, myRead.paired);
					++ numTails_;
					// rSeq = rSeq.substr(0, rBegin) ;
					// myRead.read = rSeq;
				} else { // rm read
					mylog = trimlog_t (myRead.ID, 0, rLen - 1, myRead.paired);
					++ numFull_;
					// myRead.read = std::string(1, 'x');
				}
				// log = log_t(myRead.header, rBegin, rEnd, it->dir);
				return true;
			} else { // internal matching
				if (rEnd - rBegin + 1 >= minInternalMSz_) { // rm read
					//myRead.read = std::string(1, 'x');
					//log = log_t(myRead.header, rBegin, rEnd, it->dir);
					mylog = trimlog_t (myRead.ID, 0, rLen - 1, myRead.paired);
					++ numFull_;
					return true;
				}
			}
		} // if ... else
	}
	return false;
}


/*
 * output log file and cout some statistics
 */
void Trimmer::output () {

	// write short logs to file
	std::string logFileNm = _oDir + "trim.log";

	std::ofstream olHandle(logFileNm.c_str());
	if (!olHandle.good()) {
		warning ("\tCan't create trim short log file: " + logFileNm);
	} else {
		std::cout << "\tWriting trim logs to file: " << logFileNm << "\n";
		write2log (olHandle, _logs);
		olHandle.close();
	}
	// print out
	std::cout << "\tNumber of Reads Trimmed by Heads, Tails, Full, Total:\n";
	int total = numHeads_ + numTails_ + numFull_;
	std::cout << "\t" << numHeads_ << "(" << 100.0 * (numHeads_) /total
			<< "%total)\t"
		<< numTails_	<< "(" << 100.0 * numTails_/total << "%total)\t"
		<< numFull_ << "("
		<< 100.0 * (numFull_)/total << "%total)\t" << total << "\n";
	std::cout << "\tNumber of low complexity reads: " << numlowComplexity_ << "\n";

	int pairTrimCnt = 0;
	int idx = 1;
	while (idx < (int) _logs.size()){
		int curID = _logs[idx].readID, preID = _logs[idx-1].readID;
		if (preID == (curID - 1)) {
			if (_logs[idx].paired && _logs[idx].readID % 2 == 1 &&
				_logs[idx - 1].paired && _logs[idx - 1].readID % 2 == 0) {
				++ pairTrimCnt;
				idx += 2;
			}
			else ++ idx;
		} else ++ idx;
	}
	std::cout << "\tNumber of trimmed read-pairs: " << pairTrimCnt
			<< " ( " << 100.0 * pairTrimCnt * 2/ total << " %total)\n";

} // output

void Trimmer::write2log (std::ofstream& olHandle,
		const std::vector<trimlog_t>& trimLogs) {

	for (int i = 0; i < (int) trimLogs.size(); ++ i) {
		olHandle << trimLogs[i].readID << "\t"
				 << (int) trimLogs[i].from << "\t"
				 << (int) trimLogs[i].to << "\n";
	}
} // write2log
