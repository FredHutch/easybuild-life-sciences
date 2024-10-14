//============================================================================
// Project     : DivAnalysis
// Name        : Trimmer.h
// Author      : Xiao Yang
// Created on  : Aug 2, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef TRIMMER_H_
#define TRIMMER_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include "readEntry.h"
#include "jaz/string_add.hpp"

typedef struct VectorSeedInfo {
	VectorSeedInfo (int vId, int pos, bool d):
		vecID(vId), start (pos), dir (d){}
	VectorSeedInfo (){}
	int vecID;
	int start;
	bool dir;
} seed_t;

/* the log of trimming, to be retained in memory
 * NOTE: [from] and [to] are NOT the matching position for vectors,
 * they are adjusted w.r.t. overhang, min Read length etc.
 */
struct trimlog_t {
	trimlog_t(){ readID = -1; }
	trimlog_t (int rID, int f, int t):
		readID(rID), from(f), to (t) {}
	trimlog_t (int rID): readID(rID) {}
	int readID;
	int from; // w.r.t. fwd strand trimming start
	int to;   // trimming end
} ;

inline bool operator<(const trimlog_t& lhs, const trimlog_t& rhs) {
    return lhs.readID < rhs.readID;
} // operator<


class Trimmer {
public:
	//Trimmer();
	Trimmer(const std::string& fileNm) {
		read_trim_log(fileNm);
	}

	/***********************************************************************
	 * @brief	Read trim logs from an input file in tab delimited format:
	 * 			[rID] [start_pos] [end_pos]
	 ***********************************************************************/
	void read_trim_log (const std::string& logNm) {
		std::ifstream iHandle (logNm.c_str());
		if (!iHandle.good()) {
			std::cout << "Can't open trim log file: " + logNm << "\n";
			exit(1);
		} else {
			std::cout << "Reading trim logs from file: " << logNm << "\n";
			std::string cur_line;
			while (std::getline(iHandle, cur_line)) {
				std::vector<std::string> elems;
				if (cur_line.empty()) continue; // skip empty line
				jaz::split ('\t', cur_line, std::back_inserter(elems));
				int rID = std::atoi(elems[0].c_str()),
					from = std::atoi(elems[1].c_str()),
					to = std::atoi(elems[2].c_str());
				_logs.push_back(trimlog_t(rID, from, to));
			}
			std::cout << "\tretrieved " << _logs.size() << " logs\n\n";
		}
	} // read_trim_log

	/***********************************************************************
	 * @brief	Identify readID in trim log [_logs], which is used to apply
	 * 			trimming to the read when applicable, the trimmed read is
	 * 			stored as [myRead]
	 **********************************************************************/
	void trim (std::string& myRead, const rentry_t& rEntry) const {

		myRead = rEntry.read;

		typedef 	std::vector<trimlog_t>::const_iterator iter_t;

		std::pair<iter_t, iter_t> bounds = std::equal_range(
				_logs.begin(), _logs.end(), trimlog_t(rEntry.ID));

		if (bounds.first != bounds.second) { // found
			if (bounds.first->from == 0) {
				// rm read
				if (bounds.first->to == (int) myRead.length() - 1) {
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

	std::vector<trimlog_t>& getLogs() { return _logs; }

private:
	std::vector<trimlog_t> _logs;
};

#endif /* TRIMMER_H_ */
