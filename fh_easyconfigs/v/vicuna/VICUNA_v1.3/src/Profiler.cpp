//============================================================================
// Project     : Diversifier
// Name        : Profiler.cpp
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#include "Profiler.h"
#include "jaz/fasta_file.hpp"


Profiler::~Profiler() {
	// TODO Auto-generated destructor stub
}

void Profiler::get_binned_rIDs (ivec_t& rIDs) {
	for (int i = 0 ; i < (int) binship_.size(); ++ i)
		rIDs.push_back(binship_[i].rID);
}

void Profiler::run (const ReadIndexer& rIndexer, const Trimmer& myTrimmer){

	if (iMSAFileName_.empty()) {
		std::cout << "No MSA file specified, Profiler Skipped\n\n";
		return;
	}

	std::cout << "Profiler ...\n\n";
	strvec_t genomes;
	readMSA (genomes);

	std::cout << "\tProfiling ... ...\n";
	profiling (genomes);

	std::cout << "\tReplicate profiles ... ...\n";
	profileDupl ();

	std::cout << "\tRead assignment ... ...\n";
	readAssignment (rIndexer, myTrimmer);

}


void Profiler::readMSA (strvec_t& seqs){

	readFasta (std::back_inserter(seqs), iMSAFileName_);

	// sanity check whether all sequences have the same length
	int length;
	if (seqs.size()) length = seqs[0].length();
	else warning ("no sequence in file "+ iMSAFileName_);

	for (int i = 1; i < (int) seqs.size(); ++ i) {
		if (length != (int) seqs[i].length()) {
			warning ("Sequences in MSA file " + iMSAFileName_ +
					" have different length");
		}
	}
}

void Profiler::profiling(const strvec_t& seqs) {

	// 2. identify boundaries of each bin for each input genome, stored
	//     in [coord]
	std::vector<std::vector<ipair_t> > coord; //[ )
	identifyBoundary (coord, seqs);

	// print out coordinates
	/*
	std::cout << "-----------------------------------------------------------\n";
	std::cout << "Coordinates of bins for each genome:\n";
	std::cout << "-----------------------------------------------------------\n";
	print_2dvecPairs(coord);
	std::cout << "-----------------------------------------------------------\n";
	std::cout << "kmer counts for each bin\n";
	std::cout << "-----------------------------------------------------------\n";
	*/

	// 2. extract k-spectrum for each bin and store in [_kmers]
	for (uint8_t bIdx = 0; bIdx < binNum_; ++ bIdx) {
		// a) extract k-spectrum for bIdx's bin
		std::vector<kmer_t> kmers;
		for (int gIdx = 0; gIdx < (int) coord.size(); ++ gIdx) {
			std::string seq = seqs[gIdx].substr(coord[gIdx][bIdx].first,
					coord[gIdx][bIdx].second - coord[gIdx][bIdx].first + 1);
			xny::remove_gap(seq);

			xny::get_bitkmer <std::back_insert_iterator<std::vector<kmer_t> >,
							kmer_t> (seq, std::back_inserter(kmers), K_, 0);
			// KmerExtr (seq, kmers, k, 0);
		}
		std::sort(kmers.begin(), kmers.end());

		// b) remove kmers with low frequency (= 1 currently)

		std::vector<kmer_t> FreqKmers;
		if (seqs.size() < 3) { // if there are only a few input genomes
			FreqKmers = kmers;
		}
		else {
			bool flag = false;
			for (int i = 1; i < (int) kmers.size(); ++ i) {
				if (kmers[i] == kmers[i - 1]) {
					flag = true;
				} else {
					if (flag) {
						FreqKmers.push_back(kmers[i - 1]);
						flag = false;
					}
				}
			}
			if (flag) FreqKmers.push_back(kmers.back());
		}
		// c) store in _kmers
		for (int i = 0; i < (int) FreqKmers.size(); ++ i) {
			_kmers.push_back(std::pair<kmer_t, uint8_t>(FreqKmers[i], bIdx));
			// stores reverse complimentary as well in order to speedup mapping process
			_kmers.push_back(std::pair<kmer_t, uint8_t>(
				xny::get_rvc_bits<kmer_t> (FreqKmers[i], K_), bIdx));
		}
	}

	std::sort (_kmers.begin(), _kmers.end(), compKP<kmer_t>());
	// print out
	std::cout << "\t_kmers size = " << _kmers.size() << "\n";

	// generate histogram of # of bins each kmer belongs to, in the form of
	// #kmers #bins (each of the kmer belongs to)
	std::map<kmer_t, int> kCnt;
	for (int i = 0; i < (int) _kmers.size(); ++ i) {
		std::map<kmer_t, int>::iterator it (kCnt.find(_kmers[i].first));
		if (it != kCnt.end()) {
			it->second ++;
		} else kCnt[_kmers[i].first] = 1;
	}
	ivec_t kbHisto;
	for (std::map<kmer_t, int>::iterator it = kCnt.begin();
			it != kCnt.end(); ++ it) {
		kbHisto.push_back(it->second);
	}
	std::cout << "\t-------------------------------------------------\n";
	std::cout << "\tHistogram for #bins each kmer belongs to\n";
	std::cout <<"\t\t#bins\t#kmer\n";
	histogram (kbHisto, &std::cout);
	std::cout << "\t-------------------------------------------------\n";

	//3. [TODO?] clean _kprofile: remove a kmer if it belongs to too many bins
}

// kmers that span boundaries of bins are assigned to the next bin
std::vector<ipair_t> Profiler::adjust (std::vector<ipair_t> coord,
		const std::string& seq) {

	for (int i = 1; i < (int) coord.size(); ++ i) {
		int j = coord[i].first - 1;
		int cnt = 0;
		for (; j > coord[i-1].first; -- j) {
			if (seq.at(j) != '-') ++ cnt;
			if (cnt == K_ - 1) {
				coord[i].first = j;
				break;
			}
		}
	}
	return coord;
}

/*  1) identify the longest genome after removing gaps;
 *  2) calculate nt blockSize
 *	3) identify _initCoord according to this blockSize (disregard gaps)
 */
void Profiler::initCoord (const strvec_t& seqs) {
	std::string myRef;
	int idx, maxLen = 0;
	for (int i = 0; i < (int) seqs.size(); ++ i) {
		myRef = seqs[i];
		xny::remove_gap(myRef);
		int len = (int) myRef.length();
		if (len > maxLen) {
			maxLen = len;
			idx = i;
		}
	}

	myRef = seqs[idx];
	int refLen = (int) myRef.length();
	int blockSize = maxLen/(int) binNum_;
	int idxRef = 0, ntCounter = 0;
	int crdBegin = 0;
	while (1) {
		while (ntCounter < blockSize && idxRef < myRef.length()) {
			if (myRef.at(idxRef) != '-') ++ ntCounter;
			++ idxRef;
		}
		if (idxRef == refLen) {
			_initCoord.rbegin()->second = refLen - 1;
			break;
		}
		_initCoord.push_back(ipair_t(crdBegin, idxRef));
		crdBegin = idxRef + 1;
		ntCounter = 0;
	}
	// print out
	/*
	std::cout << "MSA length = " << seqs[0].length() << "\n";
	std::cout << "Init Coord:\n";
	for (int i = 0; i < _initCoord.size(); ++ i) {
		std::cout << _initCoord[i].first << "\t" << _initCoord[i].second <<"\n";
	}
	std::cout << "\ndone!\n";
	exit(1);
	*/
}
void Profiler::identifyBoundary (std::vector<std::vector<ipair_t> >& coord,
		const std::vector<std::string>& seqs){

	int genomeLen = 0;
	if (seqs.size()) genomeLen = seqs[0].length();

	initCoord (seqs);

	//1. initialize coordinates for each genome, [ )
	/*
	int blockSize = genomeLen / binNum;

	for (int i = 0; i < binNum; ++ i) {
		_initCoord.push_back(ipair_t(i * blockSize, (i + 1) * blockSize - 1));
	}
	// set the second point of the last bin to the end of the genome
	if (_initCoord.back().second >= genomeLen) {
		abording ("BinNum is too large, exiting ...");
	} else {
		_initCoord.rbegin()->second = (genomeLen - 1);
	}
	*/
	// 2. adjust for each genome
	for (int i = 0; i < (int) seqs.size(); ++ i) {
		coord.push_back( adjust (_initCoord, seqs[i]) );
	}
}


void Profiler::profileDupl (){
	masking (_masks, blockNum_, K_, HD_);
	for (int i = 0; i < (int) _masks.size(); ++ i) {
		std::sort (_kmers.begin(), _kmers.end(), compKT<kmer_t>(_masks[i]));
		_duplTb.insert(_duplTb.end(), _kmers.begin(), _kmers.end());
	}
	_kmers.clear(); // release memory

	// print out
	std::cout << "\t\tTotal entry of _duplTb: " << (int) _duplTb.size() << "\n";
}

/* resolve ambiguities in paired-read mapping using library Size information
 have to assign read uniquely!! give higher priority to paired read info
 Consider each paired read
  1) If multiple choices satisfy libSize constraint, pick the pair with the
    highest average weight
  2) Otherwise, pick the Bin with highest weight for individual read
 */
void Profiler::rslvAmbig (int libSizeub, const ReadIndexer& rIndexer) {

	ivec_t binStat (_initCoord.size(), 0); // binStat[i]: # of reads assigned to bin i

	// consider binship of read pairs
	// for (int i = 0; i < (int) _binship.size(); ++ i) {
	for (int i = 0; i < (int) binship_.size(); ++ i) {

		uint8_t bin1, bin2, avgW = 0, maxW1 = 0, maxW2 = 0;

		bool isPaired = false; // indicate whether it is a paired read mapping

		if ((i + 1) < (int) binship_.size() &&
			rIndexer.isPaired (binship_[i].rID) &&
			(binship_[i+1].rID - binship_[i].rID) == 1 &&
			(binship_[i].rID % 2) == 0) { // two adjacent reads are paired

			for (int r1 = 0; r1 < (int) binship_[i].vec_bin_weight.size();
					r1 ++) {
				uint8_t b1 = binship_[i].vec_bin_weight[r1].first;
				uint8_t w1 = binship_[i].vec_bin_weight[r1].second;
				if (w1 > maxW1) {
					maxW1 = w1;
					bin1 = b1;
				}
				for (int r2 = 0; r2 < (int)
					binship_[i + 1].vec_bin_weight.size(); r2 ++) {
					uint8_t b2 = binship_[i + 1].vec_bin_weight[r2].first;
					uint8_t w2 = binship_[i + 1].vec_bin_weight[r2].second;
					if (w2 > maxW2) {
						maxW2 = w2;
						bin2 = b2;
					}
					if ( (b1 <= b2 && (_initCoord[b2].first - _initCoord[b1].second
								< libSizeub) ) ||
						(b1 > b2 &&	(_initCoord[b1].first - _initCoord[b2].second
								< libSizeub) ) ) {
						isPaired = true;
						uint8_t w = (w2 + w1)/2;
						if (w > avgW){
							bin1 = b1;
							bin2 = b2;
							avgW = w;
							maxW1 = maxW2 = 101;
						}
					}// if
				}// for r2
			}// for r1

			binship_[i].vec_bin_weight =
					std::vector<u8pair_t> (1, u8pair_t (bin1, isPaired));
			++ binStat[bin1];

			binship_[i + 1].vec_bin_weight =
					std::vector <u8pair_t> (1, u8pair_t(bin2, isPaired));
			++ binStat[bin2];

			++ i;

		} else { // current read is not paired with the next one
			for (int r1 = 0; r1 < (int) binship_[i].vec_bin_weight.size();
					r1 ++) {
				uint8_t b1 = binship_[i].vec_bin_weight[r1].first;
				uint8_t w1 = binship_[i].vec_bin_weight[r1].second;
				if (w1 > maxW1) {
					maxW1 = w1;
					bin1 = b1;
				}
			}
			binship_[i].vec_bin_weight =
					std::vector <u8pair_t> (1, u8pair_t(bin1, isPaired));
			++ binStat[bin1];
		}
	} // for i

	// print out binStat
	std::cout << "\t\t# of reads assigned to each bin: \n";
	std::cout << "\t\tBin\t" << "#reads\n";
	for (int i = 0; i < (int) binStat.size(); ++ i) {
		std::cout << "\t\t" << i << "\t" << binStat[i] << "\n";
	}
	std::cout << "\t\t---------------------------------------------\n";
} // rslvAmbig

void Profiler::mapping (const std::vector<read_t>& reads,
						const Trimmer& myTrimmer) {

	int preSize = reads.size();
	std::vector<bin_t> localBinship (preSize, bin_t());

	#pragma omp parallel for shared (localBinship)

	for (int i = 0; i < preSize; ++ i) {

		//0. Read Trimming if applicable (including low complexity reads)

		std::string myRead;
		int rID = reads[i].ID;

		myTrimmer.trim(myRead, reads[i]);

		if (myRead.empty()) {
			localBinship[i] = bin_t();
			continue;
		}

		//1.  for each read, map each kmer to bins, store in [rpos2bins]

		std::vector<kmer_t> rkmers;
		xny::get_bitkmer <std::back_insert_iterator<std::vector<kmer_t> >,
				kmer_t> (myRead, std::back_inserter(rkmers), K_, 0);

		// search each kmer & rvc of the read to identify their binship
		std::vector<std::set<uint8_t> > rpos2bins (rkmers.size(),
				std::set<uint8_t>());
		typedef std::vector<std::pair<kmer_t, uint8_t> >::iterator vuit_t;
		vuit_t IterBgn = _duplTb.begin();

		for (int j = 0; j < (int) rkmers.size(); ++ j) { // for each kmer

			kmer_t ID = rkmers[j]; // no need to search for rvc since
			                      //_duplTb stores both fwd and rvc strands
			uint32_t unitTbSize = _duplTb.size()/_masks.size();

			std::pair<vuit_t, vuit_t> bounds;

			for (int m = 0; m < (int) _masks.size(); ++ m) { // each mask

				bounds = std::equal_range(IterBgn + m * unitTbSize,
						IterBgn + (m + 1) * unitTbSize,
						std::pair<kmer_t, uint8_t> (ID, 0),
						compKT<kmer_t>(_masks[m]));
				for (vuit_t it = bounds.first; it != bounds.second; ++ it) {
					// identify Hamming distance
					kmer_t tmp = (*it).first ^ ID;
	                int hd = 0;
	                while (tmp) {
	                    if (tmp & 0x3) hd ++;
	                    if (hd > HD_) {
	                    		hd = HD_ + 1;
	                    		break;
	                    }
	                    tmp >>= 2;
	                }
	                if (hd <= HD_) rpos2bins[j].insert(it->second);

				} // for it
			} // for m (each mask)
		} // for j (each kmer)

		// 2. reverse mapping, [bin2pos]: bin -> pos1, pos2...
		iivec_t bin2pos (_initCoord.size(), ivec_t());
		for (int j = 0; j < (int) rpos2bins.size(); ++ j) {
			std::set<uint8_t>::iterator it (rpos2bins[j].begin());
			for (; it != rpos2bins[j].end(); ++ it){
				bin2pos[(*it)].push_back(j);
			}
		}

		// 3. calculate the cover of each bin for read i and retain
		// bins that cover at least 75% of the read

		//std::vector<ipair_t> candidate;
		std::vector<u8pair_t> candidate;

		float perc_cov = minSpan_ * myRead.length() / 100.0;

		int preCov = 0; // coverage of the previous adjacent bin
		for (uint8_t j = 0; j < bin2pos.size(); ++ j) {
			if (bin2pos[j].size()) {
				// calculate coverage for bin j
				int cov = K_;
				for (int l = 1; l < (int) bin2pos[j].size(); ++ l) {
					int dist = bin2pos[j][l] - bin2pos[j][l-1];
					if (dist < K_) cov += dist;
					else cov += K_;
				}
				// either map directly or check for boundary reads
				if (cov >= perc_cov) {
					candidate.push_back(u8pair_t (j, 100 * cov/myRead.length()));
				} else { // if cov criterion fail, check for boundary reads
					if (preCov < perc_cov && (preCov + cov) >= myRead.length()) {
						if (preCov > cov) {
							candidate.push_back(u8pair_t (j-1,
								100 * preCov/myRead.length()));
						}
						else {
							candidate.push_back(u8pair_t (j,
								100 * cov/myRead.length()));
						}
					}
				}

				preCov = cov;

			} else {
				preCov = 0;
			}
		}
		// push back no matter whether candidate is empty or not
		localBinship[i] = bin_t(rID, candidate);

		// print out
		/*
		std::cout << "\tbin\tcover\n";
		for (int j = 0; j < binCover.size(); ++ j){
			if (binCover[j] > 0)
			std::cout << "\t" << j << ":" << binCover[j] << "\n";
		}
		std::cout << "---------------------------------------------------------\n";
		*/
	}

	//_binship.insert(_binship.end(), localBinship.begin(), localBinship.end());
	for (int i = 0; i < preSize; ++ i) {
		if (localBinship[i].vec_bin_weight.size()) binship_.push_back(localBinship[i]);
	}
} // mapping

/* assign paired-reads (in FileName1 and FileName2, respectively) to
 * bins using (k, d) measure
 */
void Profiler::readAssignment (const ReadIndexer& rIndexer,
							  const Trimmer& myTrimmer){

	std::vector<read_t> reads;
	int counter = 0;
	ivec_t rIDs;
	for (int i = 0; i < (int) rIndexer.size(); ++ i) {
		rIDs.push_back(i);
		counter ++;
		if (counter >= batchSize_) {
			rIndexer.getReads(std::back_inserter(reads), rIDs.begin(),
					rIDs.end());
			// process this batch
			mapping (reads, myTrimmer);
			counter = 0;
			reads.clear();
			rIDs.clear();
		}
	}
	if (counter > 0) {
		rIndexer.getReads(std::back_inserter(reads), rIDs.begin(),
				rIDs.end());
		// process the remaining
		mapping (reads, myTrimmer);
		rIDs.clear();
		reads.clear();
	}

	// clear memory
	_duplTb.clear();

	std::cout << "\t\trslv ambiguities ...\n";
	rslvAmbig (libSzUb_, rIndexer);

	outputBinnedRIDs (oRMapFileName_);
} //readAssignment

/* Store the IDs of reads that have been assigned to Bins to a file
 *
 */
void Profiler::outputBinnedRIDs(const std::string& oFileNm) {

	int total_mapped_reads = 0, paired_reads = 0;
	std::ofstream rIDMapfHandle, rMapfHandle;
	std::string oFileNm2;
	if (!oFileNm.empty()) {
		oFileNm2 = oFileNm + ".record.txt";
    		rIDMapfHandle.open(oFileNm.c_str());
    		rMapfHandle.open (oFileNm2.c_str());
    }
	if (rIDMapfHandle.good()) {
		rIDMapfHandle << "ReadID\tBin\tisPaired?\n";
	}
	int numEntry = 0;
	for (int i = 0 ; i < (int) binship_.size(); ++ i) {
		if (binship_[i].vec_bin_weight.front().second) {
			paired_reads ++;
		}
		total_mapped_reads ++;
		if (rIDMapfHandle.good() && rMapfHandle.good()) {

			rIDMapfHandle << binship_[i].rID << "\t"
					<< (int) binship_[i].vec_bin_weight.front().first
					<< "\t" << (int) binship_[i].vec_bin_weight.front().second
					<<"\n";

			rMapfHandle << binship_[i].rID << "\t";
			if (numEntry != 0 && numEntry % 10 == 0) {
				rMapfHandle << "\n";
				numEntry = 0;
			}
			numEntry ++;
		}
	}
	rMapfHandle.close();
	rIDMapfHandle.close();
	std::cout << "\t\tStat ...\n";
	std::cout << "\t\t# reads that are mapped as pairs: " << paired_reads << "\n";
	std::cout << "\t\t# mapped reads: " << total_mapped_reads << "\n";
	std::cout << "\t\t---------------------------------------------\n";
} // outputBinnedReads
