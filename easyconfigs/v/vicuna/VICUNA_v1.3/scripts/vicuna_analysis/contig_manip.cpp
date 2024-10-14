//============================================================================
// Project     : DivAnalysis
// Name        : contig_manip.cpp
// Author      : Xiao Yang
// Created on  : Dec 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================
#include <sstream>
#include "contig_manip.h"

/***********************************************************************
 * @brief	Read input contig alignment information from file [fileNm]
 **********************************************************************/
void ContigManip::readinput() {

	std::ifstream iHandle (_aln_fNm.c_str());
	if (!iHandle.good()) {
		std::cout << "Can't open contig alignment file: " + _aln_fNm << "\n";
		exit(1);
	} else {
		std::cout << "Reading contig alignment file: " << _aln_fNm << "\n";

		std::string cur_line;
		while (std::getline(iHandle, cur_line)) {
			std::vector<std::string> elems;
			if (!cur_line.empty()) {
				if (cur_line.at(0) == '>') {
					jaz::split('\t', cur_line,
							std::back_inserter(elems));
					elems[0].erase(0, 1); // remove '>'

					contig_aln_t contig;

					contig.cID = std::atoi(elems[0].c_str());
					contig.cLen = std::atoi(elems[1].c_str());

					int 	rNum = std::atoi(elems[2].c_str());

					/* read entry */
					cread_t cread;
					for (int i = 0; i < rNum; ++ i) {
						iHandle >> cread.rID >> cread.rLen
							>> cread.rName >> cread.dist_to_beg
							>> cread.is_paired >> cread.is_fwd_strand;
						int num_indel;
						iHandle >> num_indel;
						cread.indels.resize(num_indel);
						/* indels */
						for (int j = 0; j < num_indel; ++ j) {
							iHandle >> cread.indels[j].first
									>> cread.indels[j].second;
						}

						contig.rElems.push_back(cread);

					} // for (int i =

					_cdb.push_back(contig);
				}
			} // if
		} // while

		/* print out contig information */
		for (int i = 0; i < (int) _cdb.size(); ++ i) {
			std::cout << "\t" << "cID: " <<  _cdb[i].cID
					<< "\tcLen: " << _cdb[i].cLen
					<< "\trNum:" << _cdb[i].rElems.size() << "\n";
		}

	}
} // readinput

/**********************************************************************
 * @brief	Identify low-frequency-variant free consensus for
 * 			each contig
 **********************************************************************/
void ContigManip::output_contigs (const ReadIndexer& rIndexer,
		const Trimmer& myTrimmer) {

	std::cout << "output contigs \n";

	int cNum = _cdb.size();
	std::string lfv = _oDir + "/contig.lfv.fasta",
				lfv_free = _oDir + "/contig.fa",
				lfv_info = _oDir + "/lfv.info.txt";
	std::ofstream oH_lfv (lfv.c_str()),
				  oH_lfv_free (lfv_free.c_str()),
				  oH_lfv_info (lfv_info.c_str());
	if (!oH_lfv.good()) {
		std::cout << "ContigManip::id_lfv_free_contigs - err opening "
				  << lfv << ", no output will be generated\n";
	    return;
	}
	if (!oH_lfv_free.good()) {
		std::cout << "ContigManip::id_lfv_free_contigs - err opening "
				  << lfv_free << ", no output will be generated\n";
	    return;
	}
	if (!oH_lfv_info.good()) {
			std::cout << "ContigManip::id_lfv_free_contigs - err opening "
					  << lfv_info << ", no output will be generated\n";
	} else {
		oH_lfv_info << "Contig_ID\tStart\tEnd\n\n";
	}

	/* output contigs as fasta format */
	for (int i = 0; i < cNum; ++ i) {
		/* get all rIDs involved */
		ivec_t rIDs;
		get_rIDs (rIDs, i);

		/* get all reads involved and apply trimming */
		std::vector<rentry_t> reads;
		rIndexer.get_trimmed_reads(std::back_inserter(reads), rIDs.begin(),
									rIDs.end(),	myTrimmer);

		/* identify profile + consensus of contig [i] */
		std::string consensus;
		ivec_t profile;
		generate_profile (consensus, profile, i, reads);

		/* output contigs with lfvs */
		oH_lfv << ">contig_" << i << "_lfv\t" << consensus.length()
				<< " bp\n" << consensus << "\n";

		/* identify lfv of contig [i]*/
		std::vector<ipair_t> lfv;
		id_lfv (lfv, profile, i, consensus.length());
		std::string lfv_free_contig;
		get_lfv_free_contig (lfv_free_contig, consensus, lfv);

		/* output lfv information */
		for (int j = 0; j < (int) lfv.size(); ++ j) {
			oH_lfv_info << i << "\t" << lfv[j].first
					<< "\t" << lfv[j].second << "\n";
		}

		/* output contigs w/o lfvs */
		oH_lfv_free << ">contig_" << i << "\t" << lfv_free_contig.length()
				<< " bp\n" << lfv_free_contig << "\n";
	}

	oH_lfv.close();
	oH_lfv_free.close();
	oH_lfv_info.close();
} // output_contigs

/* identify all lfv(s) */
void ContigManip::id_lfv (std::vector<ipair_t>& lfv,
				const ivec_t& profile, int cID, int cLen) {

	// position and type: drop or rise (left bracket or right bracket)
	ivec_t weights;
	std::vector<std::pair<int, char> > anchors;
	int pre_w = 0;
	for (int i = 0; i < 4; ++ i) pre_w += profile[i];
	pre_w = std::max(1, pre_w); // prevent dividing 0
	weights.push_back(pre_w);
	for (int col = 1; col < cLen; ++ col) {
		int cur_w = 0;
		for (int i = 0; i < 4; ++ i) cur_w += profile[4 * col + i];
		cur_w = std::max(1, cur_w);
		weights.push_back(cur_w);
		if (100 * cur_w / pre_w <= _lfv_freq) {
			// coverage dropping
			anchors.push_back(std::pair<int, char> (col, 'l'));
		} else if (100 * pre_w / cur_w <= _lfv_freq) {
			anchors.push_back(std::pair<int, char> (col - 1, 'r'));
		}
		pre_w = cur_w;
	}

	int index = (int) anchors.size() - 1;
	while (index >= 0) {
		if (anchors[index].second == 'r') {
			bool paired = false;
			int nested_r = 0;

			for (int i = index - 1; i >= 0; -- i) {
				int l = anchors[index].first - anchors[i].first + 1;
				if (l > _lfv_max_len) break;
				if (anchors[i].second == 'l') {
					if (nested_r != 0) -- nested_r;
					else {
						int idx_e = anchors[index].first + 1;
						int idx_s = anchors[i].first - 1;
						int r_weight = weights[idx_e];
						int l_weight = weights[idx_s];
						bool low_frq = true;
						for (int j = idx_s + 2; j <= idx_e - 2; ++ j) {
							if (100 * weights[j]/l_weight > _lfv_freq
							|| 100 * weights[j]/r_weight > _lfv_freq){
								low_frq = false;
								break;
							}
						}
						if (! low_frq) break;
						lfv.push_back(ipair_t(anchors[i].first, l));
						index = i - 1;
						paired = true;
						break;
					}
				} else { // nested 'r'
					++ nested_r;
				}
			}
			if (!paired) -- index;
		} else --index;
	}

	if (lfv.size()) std::reverse(lfv.begin(), lfv.end());
}

/* @brief	Get the consensus sequence [rc] after
 * 			removing low frq variants
 */
void ContigManip::get_lfv_free_contig (std::string& rcons,
		const std::string& consensus, const std::vector<ipair_t>& lfv) {

	int c_len = consensus.length();
	rcons = "";
	if (lfv.size()) {
		int start = 0;
		for (int i = 0; i < (int) lfv.size(); ++ i) {
			rcons += consensus.substr(start, lfv[i].first - start);
			start = lfv[i].first + lfv[i].second;
		}
		rcons += consensus.substr(start, c_len - start);
	} else rcons = consensus;

} // get_lfv_free_contig


void ContigManip::output (const ReadIndexer& rIndexer,
		const Trimmer& myTrimmer) {

	std::stringstream ss;

	/* output all contig alignments */
	for (int i = 0; i < (int) _cdb.size(); ++ i) {
		ss.str("");
		ss << i;
		std::string fileNm = ss.str();
		fileNm = _oDir + "/contig." + fileNm + ".txt";
		std::ofstream oHandle (fileNm.c_str());
		if (!oHandle.good()) {
			std::cout << "ContigManip::output -- err opening "
					  << fileNm << ", no output will be generated\n";
		    return;
		} else {
			std::cout << "\noutput to file: " << fileNm << "\n";
			print_aln (oHandle, i, rIndexer, myTrimmer, true, true);
			oHandle.close();
		}
	}

	/* print out selected regions */
	for (int i = 0; i < _num_region; ++ i) {
		int cID =  _regions[i].cID;
		ss.str("");
		ss << cID;
		std::string fileNm = ss.str();
		fileNm = _oDir + "/contig." + fileNm + ".";
		ss.str("");
		ss << _regions[i].start;
		fileNm += ss.str() + ".";
		ss.str("");
		ss << _regions[i].end;
		fileNm += ss.str() + ".txt";
		std::ofstream oHandle (fileNm.c_str());
		if (!oHandle.good()) {
			std::cout << "ContigManip::output -- err opening "
					  << fileNm << ", no output will be generated\n";
		    return;
		} else {
			if (cID < 0 || cID >= (int) _cdb.size()) {
				std::cout <<
					"[WARNING]: contig[" << cID << "] does not exist!\n";
			} else {
				std::cout << "\noutput to file: " << fileNm << "\n";

				print_aln (oHandle, cID, rIndexer, myTrimmer, true, true,
						_regions[i].start, _regions[i].end);
			}
		}

	} // for (int i

} // output

/***********************************************************************
 * @brief	Print read layout wrt contig [cID] to file [oHandle].
 * 		    If [start] && [end] != -1, then print out read layout
 * 		    between [start] and [end] positions only
 **********************************************************************/
void ContigManip::print_aln(std::ofstream& oHandle,
		int cID, const ReadIndexer& rIndexer, const Trimmer& myTrimmer,
		bool to_print_profile, bool to_print_reads, int start, int end) {

	/* get all rIDs involved */
	ivec_t rIDs;
	get_rIDs (rIDs, cID);

	/* get all reads involved and apply trimming */
	std::vector<rentry_t> reads;
	rIndexer.get_trimmed_reads(std::back_inserter(reads), rIDs.begin(),
								rIDs.end(),	myTrimmer);

	/* identify profile + consensus of contig */
	std::string consensus;
	ivec_t profile;
	generate_profile (consensus, profile, cID, reads);

	/* */
	if (to_print_profile) print_profile (oHandle, consensus, profile);
	if (to_print_reads) {
		if (start == -1 && end == -1) {
			if (reads.size() >= 50000) {
				std::cout << "\trNum = " << reads.size() << " is too large,"
						<< "\t reads are not printed\n";
				return;
			}
		} else print_reads (oHandle, cID, reads, start, end);
	}

} // print_aln

void ContigManip::generate_profile (std::string& consensus,
	ivec_t& profile, int cID, const std::vector<rentry_t>& reads) {

	int cLen = _cdb[cID].cLen;
	consensus = std::string(cLen, 'N');
	profile = ivec_t (cLen * 4, 0); // profile_[4*col + row], row: acgt
	int rNum = _cdb[cID].rElems.size();
	for (int i = 0; i < rNum; ++ i) {
		update_profile(consensus, profile, _cdb[cID].rElems[i],
					   reads[i].read, 1);
	}
}

void ContigManip::update_profile (std::string& consensus, ivec_t& profile,
  cread_t r_aln, const std::string& r_seq, int add) {

	std::string read = r_seq;
	if (!r_aln.is_fwd_strand) xny::rvc_str(read);
	int ridx = 0;
	int col = r_aln.dist_to_beg;
	for (int i = 0; i < (int) r_aln.indels.size(); ++ i) {
		int gap_start = r_aln.indels[i].first;
		int gap_len = r_aln.indels[i].second;
		for (; ridx < gap_start; ++ ridx, ++ col) {
			if (col*4 + 3 > (int) profile.size()) {
				std::cout << "Contig::create_profile p1 SC failed\n";
				exit(1);
			}
			update_profile_column (consensus, profile, col, read.at(ridx), add);
		}
		col += gap_len;
	}
	for (; ridx < r_aln.rLen; ++ ridx, ++ col) {
		if (col*4 + 3 > (int) profile.size()) {
			std::cout << "Contig::create_profile p2 SC failed\n";
			exit(1);
		}
		update_profile_column (consensus, profile, col, read.at(ridx), add);
	}
} // Contig::update_profile

void ContigManip::update_profile_column (std::string& consensus,
		ivec_t& profile, int col, char newchar, int value)  {

	int consensus_row = xny::char2bits(consensus.at(col));

	int new_row = xny::char2bits(newchar);

	if (new_row != -1) {

		profile[col*4 + new_row] += value;

		if (consensus_row == -1) {
			consensus.at(col) = std::toupper(newchar);
			return;
		}

		if (consensus_row != new_row) {
			if (	profile[col*4 + new_row] > profile[col*4 + consensus_row])
				consensus.at(col) = std::toupper(newchar);
		} else {
			if (value < 0) {
				for (int r = 0; r < 4; ++ r) {
					if (profile[col*4 + r] > profile[col*4 + new_row]) {
						consensus.at(col) = std::toupper(xny::bits2char(r));
						break;
					}
				} // for
			}
		}
	} // if (row != -1)
} // update_profile_column

void ContigManip::print_profile(std::ofstream& oHandle,
		const std::string& consensus, const ivec_t& profile) {

	for (int row = 0; row < 4; ++ row) {
		for (int col = 0; col < (int) profile.size()/4; ++ col) {
			oHandle << std::setw(7) << profile[4*col + row];
		}
		oHandle << "\n";
	}
	for (int i = 0; i < (int) consensus.length(); ++ i) {
		oHandle << std::setw(7) << consensus.at(i);
	}
	oHandle << "\n";
	for (int i = 0; i < (int) consensus.length(); ++ i) {
		oHandle << std::setw(7) << i;
	}
	oHandle << "\n\n";
	oHandle << consensus << "\n\n";
} // print_profile


void ContigManip::print_reads (std::ofstream& oHandle, int cID,
		const std::vector<rentry_t>& reads, int start, int end){

	oHandle << "Region: [ " << start << ", " << end << " ]\n\n";

	int rNum = reads.size();
	for 	(int i = 0; i < rNum; ++ i) {
		cread_t r_aln = _cdb[cID].rElems[i];
		rentry_t r_entry = reads[i];

		if (start != -1 && end != -1) {
			int read_end = r_aln.dist_to_beg + r_aln.rLen +
							r_aln.indel_length();
			if (read_end < start || r_aln.dist_to_beg > end) continue;
		}

		for (int j = 0; j < r_aln.dist_to_beg; ++ j) oHandle << "-";

		std::string r_seq = r_entry.read;

		if (!r_aln.is_fwd_strand) xny::rvc_str (r_seq);

		int ridx = 0;
		for (int j = 0; j < (int) r_aln.indels.size(); ++ j) {
			int gap_start = r_aln.indels[j].first;
			int gap_len = r_aln.indels[j].second;
			oHandle << r_seq.substr(ridx, gap_start - ridx);
			for (int g = 0; g < gap_len; ++ g) oHandle << "-";
			ridx = gap_start;
		}
		oHandle << r_seq.substr(ridx, r_aln.rLen - ridx) << "\t";

		if (r_aln.is_fwd_strand) oHandle << "+\t";
		else oHandle << "-\t";
		oHandle << r_aln.rID << "\t" << r_aln.dist_to_beg << "\n";
				//<< "\t" << r_entry.header << "\n";
	}

} // print_read
