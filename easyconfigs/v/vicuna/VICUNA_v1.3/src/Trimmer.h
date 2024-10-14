//============================================================================
// Project     : Diversifier
// Name        : Trimmer.h
// Author      : Xiao Yang
// Created on  : Aug 2, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef TRIMMER_H_
#define TRIMMER_H_
#include "jaz/fasta_file.hpp"
#include "xny/seq_manip.hpp"
#include "xny/seq_pair_manip.hpp"
#include "xutil.h"
#include "Parameter.h"

class ReadIndexer;

typedef struct ReadEntry read_t;

struct seed_t {
	seed_t (int vId, int pos, bool d):
		vecID(vId), start (pos), dir (d){}
	seed_t (){}
	int vecID;
	int start;
	bool dir;
} ;

/* the log of trimming, to be retained in memory
 * NOTE: [from] and [to] are NOT the matching position for vectors,
 * they are adjusted w.r.t. overhang, min Read length etc.
 */
struct trimlog_t {
	trimlog_t(){ readID = -1; }
	trimlog_t (int rID, int f, int t, bool flag):
		readID(rID), from(f), to (t), paired (flag) {}
	trimlog_t (int rID): readID(rID) {}
	int readID;
	uint8_t from; // w.r.t. fwd strand trimming start
	uint8_t to;   // trimming end
	bool paired;
};

inline bool operator<(const trimlog_t& lhs, const trimlog_t& rhs) {
    return lhs.readID < rhs.readID;
} // operator<


class Trimmer {
public:
	//Trimmer();
	Trimmer(const Parameter& myPara);
	virtual ~Trimmer();
	void run (const ReadIndexer& rIndexer);

	// apply trimming to a read
	void trim(std::string& myRead, const read_t& read) const;
	std::vector<trimlog_t>& getLogs() { return _logs; }

private:
	std::vector<trimlog_t> _logs;
	// parameters from input
	int maxOverhang_;
	int minInternalMSz_;
	int minMSz_;
	int minReadSz_;
	int _batchsz;
	std::string _vecName;
	std::string _logName;
	std::string _oDir;

	// some statistics of trimming
	int numHeads_;
	int numTails_;
	int numFull_;
	int numlowComplexity_;

	void read_trim_log ();

	void processFile (const strvec_t& vectors, const std::string& fileNm,
			const std::string& trimFileNm);
	void vBatchTrim (const std::vector<read_t>& reads, const strvec_t& vectors);
	bool trimmable (trimlog_t& log, const std::vector<seed_t>& kinfo,
			const strvec_t& vectors, const read_t& myRead, int rPos);

	void output ();

	void write2log (std::ofstream& olHandle,
			const std::vector<trimlog_t>& trimSLogs);

	//TODO void qualityTrimming();
};

#endif /* TRIMMER_H_ */
