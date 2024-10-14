//============================================================================
// Project     : Diversifier
// Name        : Contig.cpp
// Author      : Xiao Yang
// Created on  : Aug 26, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#include "Contig.h"

/* @brief	Initialize contig using super-shingle information. Reads are
 * 			anchored by a kmer that contributes to the min_ss
 */
void Contig::init (const ivec_t& shingle_indices,
				   const std::vector<shingle_t>& shingles,
				   const std::vector<read_t>& reads){

	int numReads = shingle_indices.size();

	if (numReads < 2) {
		std::cout << numReads << "\n";
		std::cout << "Contig::init p1: sanity check failed\n";
		exit(1);
	}

	// identify the super-shingle shared among all reads
	bool flag_found = true;
	uint32_t common_ss = shingles[shingle_indices[0]].fwd_ss;
	for (int i = 1; i < numReads; ++ i) {
		if (common_ss != shingles[shingle_indices[i]].fwd_ss &&
				common_ss != shingles[shingle_indices[i]].rv_ss) {
			flag_found = false;
		}
	}
	if (!flag_found) common_ss = shingles[shingle_indices[0]].rv_ss;

	// initialize each [cread] based on [shingles] information,
	// where read sequences are retrieved from [reads]
	std::vector<cread_t> tmp_db;
	typedef std::vector<read_t>::const_iterator cit;
	int max_min_pos = 0;
	for (int i = 0; i < numReads; ++ i) {

		int idx = shingle_indices[i];

		cread_t cread;
		cread.rID = shingles[idx].rID;
		cread.is_fwd_strand =
				(shingles[idx].fwd_ss == common_ss) ? true :false;
		cread.dist_to_beg = cread.is_fwd_strand ?
				shingles[idx].fwd_min_pos: shingles[idx].rv_min_pos;
		cread.is_paired = shingles[idx].is_paired;
		cread.indels = std::vector<u8pair_t> ();
		max_min_pos = std::max (cread.dist_to_beg, max_min_pos);

		// get read using rID by binary search reads
		std::pair<cit, cit> bounds =
			std::equal_range(reads.begin(), reads.end(), read_t(cread.rID));
		if (bounds.first != bounds.second) { // found the read
			cread.read_length = bounds.first->read.length();
			tmp_db.push_back(cread);
		} else {
			std::cout << "Contig::init p2: sanity check failed\n";
			exit(1);
		}
	}

	/* update dist_to_beg by deduction using max_min_pos, and length_
	 * insert to db_ afterwards
	 */
	for (int i = 0; i < (int) tmp_db.size(); ++ i) {
		{ // debug
			if (max_min_pos < tmp_db[i].dist_to_beg) {
				std::cout << "Contig::init p2!!!!!\n";
			}
			if ( tmp_db[i].indel_length() != 0) {
				std::cout << "Contig::init p3 !!!\n";
			}
		}
		tmp_db[i].dist_to_beg = max_min_pos - tmp_db[i].dist_to_beg;
		length_ = std::max(length_,  tmp_db[i].read_length +
				 tmp_db[i].dist_to_beg +  tmp_db[i].indel_length());
		db_.insert(tmp_db[i]);
	}

} // Contig::init

/*
 * split initial contigs to multiple ones when possible,
 * and store resulting contigs in [contigs]
 * return true either in case of no splitting or no consensus can be formed
 */
bool Contig::validate (const std::set<ipair_t, comp_ipair>& rIndices,
 const std::vector<read_t>& reads, std::set<ipair_t, comp_ipair>& diverged,
 int overhang, int divergence) {

	/* Form consensus and profile for the current [db_] */
	std::string consensus;
	ivec_t profile;
	generate_profile (consensus, profile, rIndices, reads);

	/* recursively aligning constituent reads against consensus */
	while (!verify_consensus (consensus, profile, rIndices, reads, diverged,
			overhang, divergence));

	/* no splitting */
	if (diverged.size() == 0) return true;

	/* splitting happens */
	if (db_.size() - diverged.size() >= 2) return false;
	else
	{
		/* no consensus can be formed with current profile
		 * greedy approach? identify a set of reads that are similar
		 * and postponed processing remaining reads in next stage
		 */
		return true;
	}

} //Contig::validate


/*
 * map all reads back to consensus, and update consensus when
 * applicable. Store the indices of reads with large divergence
 * from consensus in [diverged_rIndices]
 *
 * 1. Penalize continuous mismatches using equation (1+2n) where,
 * 1 is the mismatch starting penalty and n is the additional #
 * of penalties (can change 2 to 1.5 or 3 etc. to decrease or
 * increase the penalized effect)
 * 2. penalize consensus columns with a single base
 */
bool Contig::verify_consensus (std::string& consensus, ivec_t& profile,
  const std::set<ipair_t, comp_ipair>& rIndices,
  const std::vector<read_t>& reads, std::set<ipair_t, comp_ipair>& diverged,
  int overhang, int divergence) {

	bool update = false;
	bool continuity = false;
	int continuity_penalty = 1;
	int diverge_threshold = divergence;
	if (db_.size() == 2) diverge_threshold /= 2;

	std::set<cread_t>::iterator it_db = db_.begin();
	for (; it_db != db_.end(); ++ it_db) {

		if (diverged.count(ipair_t(it_db->rID, 0))) continue;

		int diff = 0;

		std::string read;
		std::set<ipair_t, comp_ipair>::const_iterator
				it_s = rIndices.find(ipair_t(it_db->rID, 0));
		if (it_s != rIndices.end()) read = reads[it_s->second].read;
		else abording("Contig::verify_consensus sanity check failed!\n");


		int rlength = read.length();
		if (!it_db->is_fwd_strand) xny::rvc_str(read);

		int aligned_bases = rlength;
		uint8_t ridx = 0;
		int col = it_db->dist_to_beg;
		for (int i = 0; i < (int) it_db->indels.size(); ++ i) {

			{ // debug
				std::cout << "verify_consensus: go inside here, wrong !!!!\n";
				exit(1);
			}
			uint8_t gap_start = it_db->indels[i].first;
			uint8_t gap_len = it_db->indels[i].second;
			for (; ridx < gap_start; ++ ridx, ++ col) {
				int col_cnt = profile[4*col] + profile[4*col + 1] +
						profile[4*col + 2] +  profile[4*col + 3];
				/* penalize consensus columns with a single base */
				if (col_cnt <= 1) {
					++ diff;
					continue;
				}
				if (toupper(read.at(ridx)) != toupper(consensus.at(col))){
					if (ridx >= overhang && ridx < rlength - overhang){
						if (continuity) 	diff += (1 + continuity_penalty);
						else {
							++ diff;
							continuity = true;
						}
					}
					else {
						-- aligned_bases;
						continuity = true;
					}
				} else continuity = false;
			}
			col += gap_len;
		}
		for (; ridx < it_db->read_length; ++ ridx, ++ col) {
			int col_cnt = profile[4*col] + profile[4*col + 1] +
						profile[4*col + 2] +  profile[4*col + 3];
			/* penalize consensus columns with a single base */
			if (col_cnt <= 1) {
				++ diff;
				continue;
			}
			if (toupper(read.at(ridx)) != toupper(consensus.at(col))){
				if (ridx >= overhang && ridx < rlength - overhang){
					if (continuity) 	diff += (1 + continuity_penalty);
					else {
						++ diff;
						continuity = true;
					}
				}
				else {
					-- aligned_bases;
					continuity = true;
				}
			} else continuity = false;
		}

		if (aligned_bases <= 0) {
			std::cout << "Contig::verify_consensus sanity check failed,"
					<< "number_bases_considered = "
					<< aligned_bases << "\n";
			exit(1);
		}

		if (diff * 100.0 / aligned_bases > diverge_threshold) {
			diverged.insert(*it_s);

			update_profile (consensus, profile, it_db, read, -1);

			update = true;
		}
	}

	return update ? false : true;

} // Contig::verify_consensus


void Contig::generate_profile (std::string& consensus, ivec_t& profile,
		const std::set<ipair_t, comp_ipair>& rIndices,
		const std::vector<read_t>& reads) const{

	consensus = std::string(length_, 'N');
	profile = ivec_t (length_ * 4, 0); // profile_[4*col + row], row: acgt

	std::set<cread_t>::iterator it_db = db_.begin();
	for (; it_db != db_.end(); ++ it_db) {
		std::set<ipair_t, comp_ipair>::const_iterator
			it_s = rIndices.find (ipair_t(it_db->rID, 0));
		if (it_s != rIndices.end()) {
			update_profile(consensus, profile,
					it_db, reads[it_s->second].read, 1);
		} else abording("Contig::validate sanity check failed!\n");
	}
}

/*
 * sign > 0: addition to profile
 * sign < 0: reduction to profile
 */
void Contig::update_profile (std::string& consensus, ivec_t& profile,
  std::set<cread_t>::const_iterator it, const std::string& read_seq,
  int sign) const {
	std::string read = read_seq;
	if (!it->is_fwd_strand) xny::rvc_str(read);
	uint8_t ridx = 0;
	int col = it->dist_to_beg;
	for (int i = 0; i < (int) it->indels.size(); ++ i) {
		uint8_t gap_start = it->indels[i].first;
		uint8_t gap_len = it->indels[i].second;
		for (; ridx < gap_start; ++ ridx, ++ col) {
			if (col*4 + 3 > (int) profile.size()) {
				std::cout << "Contig::create_profile p1 SC failed\n";
				exit(1);
			}
			update_profile_column (consensus, profile, col, read.at(ridx), sign);
		}
		col += gap_len;
	}
	for (; ridx < it->read_length; ++ ridx, ++ col) {
		if (col*4 + 3 > (int) profile.size()) {
			std::cout << "Contig::create_profile p2 SC failed\n";
			exit(1);
		}
		update_profile_column (consensus, profile, col, read.at(ridx), sign);
	}
} // Contig::update_profile

/* This version takes consideration of batch processing and counting
 * consensus in a different fashion: the weight of a base is equivalent
 * to the minimum distance to either end of the read
 *
 * Procedure: looking for each rID of the contig in rIndices, if not found,
 * then this read is either processed in the previous batch or will be
 * processed in the next batch. */
void Contig::generate_weighted_profile (std::string& consensus,
	ivec_t& profile, const std::set<ipair_t, comp_ipair>& rIndices,
	const std::vector<read_t>& reads) const{

	if (!consensus.size()) { // initialize
		consensus = std::string(length_, 'N');
		profile = ivec_t (length_ * 4, 0); // profile_[4*col + row], row: acgt
	}

	std::set<cread_t>::iterator it_db = db_.begin();
	for (; it_db != db_.end(); ++ it_db) {

		std::set<ipair_t, comp_ipair>::const_iterator
			it_r = rIndices.find(ipair_t(it_db->rID, 0));

		if (it_r != rIndices.end()){ // found !

			update_weighted_profile (consensus, profile, it_db,
					reads[it_r->second].read, 1);
		}
	}
}

/* @brief	Given the existing [profile] and the corresponding [consensus],
 * 			update them once a read [read_seq] is added ([sign] = +1) or
 * 			removed. The information of this read wrt to *this contig
 * 			is specified in [it].
 */
void Contig::update_weighted_profile (std::string& consensus,
		ivec_t& profile, std::set<cread_t>::const_iterator it,
		const std::string& read_seq, int sign) const {
	std::string read = read_seq;
	if (!it->is_fwd_strand) xny::rvc_str(read);
	uint8_t ridx = 0;
	int col = it->dist_to_beg;
	for (int i = 0; i < (int) it->indels.size(); ++ i) {
		uint8_t gap_start = it->indels[i].first;
		uint8_t gap_len = it->indels[i].second;
		for (; ridx < gap_start; ++ ridx, ++ col) {
			if (col*4 + 3 > (int) profile.size()) {
				std::cout << "Contig::create_profile p1 SC failed\n";
				exit(1);
			}
		    int min_dist = std::min(ridx + 1, it->read_length - ridx);
			update_profile_column (consensus, profile, col, read.at(ridx),
					sign * min_dist);
		}
		col += gap_len;
	}
	for (; ridx < it->read_length; ++ ridx, ++ col) {
		if (col*4 + 3 > (int) profile.size()) {
			std::cout << "Contig::create_profile p2 SC failed\n";
			exit(1);
		}
		int min_dist = std::min(ridx + 1, it->read_length - ridx);
		update_profile_column (consensus, profile, col, read.at(ridx),
				sign * min_dist);
	}

} // Contig::update_weighted_profile

void Contig::update_profile_column (std::string& consensus, ivec_t& profile,
		int col, char newchar, int value) const {

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
} // Contig::update_column


/*
 *  Split current contig [rhs] into two parts, return the validated
 *  contig [lhs], and keep the remaining reads as the new [rhs]
 */
Contig Contig::split (const std::set<ipair_t, comp_ipair>& rhs_indices) {

	// get all rhsIDs used later and new distance to beg
	int lhs_min_dist_to_beg = INT_MAX, rhs_min_dist_to_beg = INT_MAX;

	std::vector<cread_t> tmp_lhs, tmp_rhs;

	for (std::set<cread_t>::iterator it = db_.begin();
			it != db_.end(); ++ it) {
		if (rhs_indices.count(ipair_t(it->rID, 0))) {
			rhs_min_dist_to_beg = std::min (rhs_min_dist_to_beg,
							it->dist_to_beg);
			tmp_rhs.push_back(*it);
		} else {
			lhs_min_dist_to_beg = std::min (lhs_min_dist_to_beg,
										it->dist_to_beg);
			tmp_lhs.push_back(*it);
		}
	}


	// update the coordinates
	update_dist_to_beg (tmp_rhs, (-1) * rhs_min_dist_to_beg);
	update_dist_to_beg (tmp_lhs, (-1) * lhs_min_dist_to_beg);

	Contig lhs;
	lhs.insert (tmp_lhs.begin(), tmp_lhs.end());
	db_.clear();
	insert (tmp_rhs.begin(), tmp_rhs.end());

	calculate_length();
	lhs.calculate_length();

	return lhs;

} //  Contig::split


/* @brief	Merge [rhs] to *this: identify each gap (pos, len) wrt *this
 * 			and rhs contig, then calculate affected reads (I for *this,
 * 			D for rhs) and make modifications.
 * @para
 * 	[algn]: specifies the alignment string [algn].seq and
 * 			the start & end positions of this alignment in both contigs
 * @assure	DDII type of alignment is not allowed
 * 			'D' or 'I' is present in the end of alignment
 */
void Contig::merge_via_overlap (const xny::align_t& algn, Contig& rhs){

	/* relative position to the front of merged contig */
	int lhs_extr_dist = 0, rhs_extr_dist = 0;
	/*
	if (algn.s0_l == 0) lhs_extr_dist = algn.s1_l;
	else rhs_extr_dist = algn.s0_l;
	*/
	if (algn.s0_l <= algn.s1_l) lhs_extr_dist = algn.s1_l - algn.s0_l;
	else rhs_extr_dist = algn.s0_l -  algn.s1_l;

	/* identify all gaps in lhs and rhs, store in vector of pairs */
	std::vector<ipair_t> lhs_gap, rhs_gap;
	int I_cnt = 0, D_cnt = 0, start;
	int I_pre = 0, D_pre = 0;
	for (int i = 0; i < (int) algn.seq.length(); ++ i) {
		switch (algn.seq.at(i)){
		case 'I':	//insert indels in lhs
			if (I_cnt == 0) start = i + algn.s0_l - I_pre;
			++ I_cnt;
			break;
		case 'D':	//insert indels in rhs
			if (D_cnt == 0) start = i + algn.s1_l - D_pre;
			++ D_cnt;
			break;
		default:
			if (I_cnt > 0) {
				lhs_gap.push_back(ipair_t(start, I_cnt));
				I_pre += I_cnt;
				I_cnt = 0;
			} else if (D_cnt > 0) {
				rhs_gap.push_back(ipair_t(start, D_cnt));
				D_pre += D_cnt;
				D_cnt = 0;
			}
			break;
		}
	} // for (int i = 0 ...

	/*
	{ // debugging
		if ((find_element(cread_t(48334)) != db_.end() && lhs_gap.size())
			|| (	rhs.find_element(cread_t(48334)) != rhs.get_db_end() &&
			rhs_gap.size())) {
			std::cout << "48334\t";
			if (lhs_gap.size()) {
				std::cout << "lhs:\n";
				for (int i = 0; i < lhs_gap.size(); ++ i) {
					std::cout << "(" << lhs_gap[i].first << ", " << lhs_gap[i].second <<")\n";
				}
			}
			if (rhs_gap.size()) {
				std::cout << "rhs:\n";
				for (int i = 0; i < rhs_gap.size(); ++ i) {
					std::cout << "(" << rhs_gap[i].first << ", " << rhs_gap[i].second <<")\n";
				}
			}
		}
		if ((find_element(cread_t(48232)) != db_.end() && lhs_gap.size())
			|| (	rhs.find_element(cread_t(48232)) != rhs.get_db_end() &&
			rhs_gap.size())) {
			std::cout << "48232\t";
			if (lhs_gap.size()) {
				std::cout << "lhs:\n";
				for (int i = 0; i < lhs_gap.size(); ++ i) {
					std::cout << "(" << lhs_gap[i].first << ", " << lhs_gap[i].second <<")\n";
				}
			}
			if (rhs_gap.size()) {
				std::cout << "rhs:\n";
				for (int i = 0; i < rhs_gap.size(); ++ i) {
					std::cout << "(" << rhs_gap[i].first << ", " << rhs_gap[i].second <<")\n";
				}
			}
		}
	}
	*/

	/* update the reads in lhs & rhs when applicable */
	update_all_reads (lhs_extr_dist, lhs_gap);
	rhs.update_all_reads (rhs_extr_dist, rhs_gap);

	/* merge */
	db_.insert(rhs.get_db_begin(), rhs.get_db_end());
	length_ = std::max(algn.s0_l, algn.s1_l) + algn.seq.length() +
	  std::max(rhs.get_length() - algn.s1_r - 1, length_ - algn.s0_r - 1);
}

/***********************************************************************
 * @brief	Set the minimum value of dist_to_beg of all reads to 0
 * 			(in case it is positive or negative), and update	contig length
 ***********************************************************************/
void Contig::reset_db () {
	std::vector<cread_t> tmp_db;
	int min_dist_to_beg = INT_MAX;
	length_ = 0;
	std::set<cread_t>::iterator it = db_.begin();
	for (; it != db_.end(); ++ it) {
		tmp_db.push_back(*it);
		min_dist_to_beg = std::min(min_dist_to_beg, it->dist_to_beg);
	}
	for (int i = 0; i < (int) tmp_db.size(); ++ i) {
		tmp_db[i].dist_to_beg -= min_dist_to_beg;
		length_ = std::max (length_,
				tmp_db[i].dist_to_beg + tmp_db[i].read_length);
	}
	db_.clear();
	db_.insert(tmp_db.begin(), tmp_db.end());

} // reset_db ()

/* during the process of merging two contigs, update the read information
 * in the selected contig db.
 * extr_dist: the distance needs to be added to all reads' dist_to_beg.
 * gaps: additional gaps required to be inserted into the current contig,
 * this information need to be incorporated into each reads' dist_to_beg
 * as well as indels.
 */

void Contig::update_all_reads (int extr_dist,
		const std::vector<ipair_t>& contig_gaps) {

	ivec_t rIDs;
	get_rIDs(rIDs);
	for (int i = 0; i < rIDs.size(); ++ i) {

		/*
		{ // debug
			if (rIDs[i] == 48334 || rIDs[i] == 48232) {
				std::cout << "here\t" << rIDs[i] << "\n";
			}
		}
		*/
		std::set<cread_t>::iterator it_db = find_element (cread_t(rIDs[i]));
		int rStart = it_db->dist_to_beg;
		int rEnd = rStart + it_db->indel_length() + it_db->read_length - 1;
		int gap_ahead = 0;
		std::vector<u8pair_t> gaps;
		for (int j = 0; j < (int) contig_gaps.size(); ++ j) {
			if (contig_gaps[j].first > rEnd){/* the read is b4 current gap */
				break;
			} else if (contig_gaps[j].first <= rStart) { /* after gap */
				gap_ahead += contig_gaps[j].second;
			} else { // gap within read
				/* the newly introduced gaps are wrt to contig positions,
				 * which need to be adjusted to individual reads, which
				 * previously may contain gaps */
				int new_gap_start = contig_gaps[j].first - rStart;
				for (int g = 0; g < (int) it_db->indels.size(); ++ g) {
					if (it_db->indels[g].first <= new_gap_start) {
						if (new_gap_start <= it_db->indels[g].first +
								it_db->indels[g].second) {
							new_gap_start = it_db->indels[g].first;
						} else {
							new_gap_start -= it_db->indels[g].second;
						}
					} else break;
				}

				gaps.push_back(u8pair_t (new_gap_start, contig_gaps[j].second));
			}
		}
		update_read (it_db, extr_dist + gap_ahead, gaps);
	}
} // update_all_reads

/* @brief	Merge [rhs] to [lhs], which share a common_rID
 * @assure	As the input [rhs].size <= [lhs].size
 */
bool Contig::merge_via_rID (Contig& rhs, int common_rID) {

	std::vector<cread_t> diff_set (rhs.number_of_reads() + db_.size(), cread_t());
	std::vector<cread_t>::iterator it = std::set_difference(rhs.get_db_begin(),
			rhs.get_db_end(), db_.begin(), db_.end(), diff_set.begin());
	diff_set.resize(int (it - diff_set.begin()));
	if (diff_set.size() == 0) {	rhs.clear(); return true;}

	/* identify which rID brings two contigs together using binary search */
	std::set<cread_t>::iterator l_iter = db_.find(cread_t(common_rID));
	std::set<cread_t>::const_iterator r_iter = rhs.find_element(cread_t(common_rID));

	if (l_iter == db_.end() || r_iter == rhs.get_db_end() ||
			l_iter->rID != common_rID || r_iter->rID != common_rID) {
		abording ("Contig::merge: no read ID found, sanity check failed\n");
	}

	bool inverse_flag =
		 (l_iter->is_fwd_strand == r_iter->is_fwd_strand) ? false: true;

	int l_dist_to_beg = l_iter->dist_to_beg,
		r_dist_to_beg = r_iter->dist_to_beg;

	if (inverse_flag) r_dist_to_beg = rhs.get_length() - r_dist_to_beg -
			r_iter->indel_length() - r_iter->read_length;

	//merge
	if (l_dist_to_beg >= r_dist_to_beg) {
		for (int i = 0; i < (int) diff_set.size(); ++ i) {
			if (inverse_flag) inverse_read (rhs.get_length(), diff_set[i]);
			diff_set[i].dist_to_beg += (l_dist_to_beg - r_dist_to_beg);
			db_.insert(diff_set[i]);
		}
	} else {

		update_dist_to_beg (r_dist_to_beg - l_dist_to_beg);
		for (int i = 0; i < (int) diff_set.size(); ++ i) {
			if (inverse_flag) inverse_read (rhs.get_length(), diff_set[i]);
			db_.insert(diff_set[i]);
		}
	}

	// new length
	int common_rlength = l_iter->read_length;
	/*
	{ // debug
		if (common_rlength > 68) {
			std::cout << "common_rlength = " << common_rlength << "\n";
			std::cout << "rID:" << l_iter->rID << "\t"
					<< "dist_to_beg" << l_iter->dist_to_beg << "\t"
					<< "#indels:" << l_iter->indels.size() << "\t"
					<< "gap length= " << (int) l_iter->indel_length() << "\n\n";
		}
	}
	*/
	length_ = std::max(l_dist_to_beg, r_dist_to_beg) + common_rlength
			 +	std::max(length_ - l_dist_to_beg - common_rlength,
			rhs.get_length() - r_dist_to_beg - common_rlength);

	rhs.clear();
	return true;

} // Contig::merge

/* update either/both dist_to_beg or/and gaps in a specific read */
void Contig::update_read (std::set<cread_t>::iterator it, int dist,
		const std::vector<u8pair_t>& gaps){
	if (dist == 0 && gaps.size() == 0) return;
	cread_t tmp = (*it);
	tmp.dist_to_beg += dist;
	if (gaps.size()) {
		tmp.indels.insert(tmp.indels.end(), gaps.begin(), gaps.end());
		/* sort gaps, then merge gaps starting at the same position*/
		std::sort(tmp.indels.begin(), tmp.indels.end(), cmp_u8pair());
		for (int i = 1; i < (int) tmp.indels.size(); ++ i) {
			if (tmp.indels[i].first == tmp.indels[i-1].first) {
				tmp.indels[i-1].second += tmp.indels[i].second;
				tmp.indels.erase(tmp.indels.begin() + i);
				-- i;
			}
		}
	}
	db_.erase(it);
	db_.insert(tmp);
}

void Contig::inverse_all_reads() {
	ivec_t rIDs;
	get_rIDs(rIDs);
	for (int i = 0; i < rIDs.size(); ++ i) {
		std::set<cread_t>::iterator it_db = find_element (cread_t(rIDs[i]));
		cread_t tmp_r = *it_db;
		db_.erase(it_db);
		inverse_read (length_, tmp_r);
		db_.insert(tmp_r);
	}
	/*
	std::set<cread_t>::iterator it_db = db_.begin();
	for (; it_db != db_.end(); ++ it_db) {
		cread_t tmp_r = *it_db;
		db_.erase(it_db);
		inverse_read (length_, tmp_r);
		db_.insert(tmp_r);
	}
	*/
} // Contig::inverse_all_reads()

/* inverse the direction of reads and update corresponding information */
void Contig::inverse_read (int contig_length, cread_t& cread) {
	cread.dist_to_beg = contig_length - cread.dist_to_beg -
			cread.indel_length() - cread.read_length;
	for (int i = 0; i < (int) cread.indels.size(); ++ i) {
		cread.indels[i].first = cread.read_length - cread.indels[i].first;
	}
	std::reverse (cread.indels.begin(), cread.indels.end());
	cread.is_fwd_strand = cread.is_fwd_strand? false: true;
}


void Contig::print_alignment (std::ofstream& oHandle,
	int contigID, const ReadIndexer& rIndexer, const Trimmer& myTrimmer) {

	ivec_t rIDs;
	get_rIDs (rIDs);
	std::vector<read_t> reads;
	rIndexer.get_trimmed_reads(std::back_inserter(reads), rIDs.begin(),
			rIDs.end(),	myTrimmer);

	std::set<ipair_t, comp_ipair> rIndices;
	for (int i = 0; i < (int) rIDs.size(); ++ i)
		rIndices.insert(ipair_t(rIDs[i], i));

	/* create profile and consensus */
	std::string consensus;
	ivec_t profile;
	generate_profile (consensus, profile, rIndices, reads);

	oHandle << "------------------------------------------------------\n";
	oHandle << "Contig alignment entry in tab delimited format:\n";
	oHandle << "header:  [>contigID][contigLen][rNum]\n";
	oHandle << "content: [rID][rLen][rName][dist_to_beg][paired?][strd]"
				"[indelNum] ([start0][len1]) ([start1][len1]) ...\n";
	oHandle << "         ...\n";
    oHandle << "------------------------------------------------------\n";

	/* consensus in fasta format */
	oHandle << ">" << contigID << "\t"
			<< consensus.length() << "\t" << rIDs.size() << "\n";
	//oHandle << consensus << "\n\n";

	/* sort reads by dist_to_beg */
	std::vector<ipair_t> dist_to_rID;
	sort_by_dist (dist_to_rID);

	for (int i = 0; i < (int) dist_to_rID.size(); ++ i) {
		int rID = dist_to_rID[i].second;
		std::set<cread_t>::iterator it_r = db_.find(cread_t(rID));
		if (it_r == db_.end())
			abording("Contig:print_reads p1 sanity check failed");

		std::vector<read_t>::const_iterator
			iter = std::lower_bound(reads.begin(), reads.end(), read_t(rID));

		if (iter == reads.end())
			abording("Contig:print_reads p2 sanity check failed");

		std::string header = iter->header;

		/* print out
		 * (rID, read_length, header, dist_to_beg, dir, num_gaps,
		 * 	[start_pos_0, len_0], [start_pos_1, len_1] ...) */
		oHandle << rID << "\t" << iter->read.length() << "\t"
				<< header << "\t" << it_r->dist_to_beg << "\t"
				<< (int) it_r->is_paired << "\t"
				<< (int) it_r->is_fwd_strand << "\t";

		int num_gaps = (int) it_r->indels.size();
		oHandle << num_gaps << "\t";

		for (int j = 0; j < num_gaps; ++ j) {
			oHandle << (int) it_r->indels[j].first << "\t" <<
					 (int) it_r->indels[j].second << "\t";
		}
		oHandle << "\n";
	}

} // print_alignment

/* @brief	To print read tiling over the contig.
 * @param
 * 	[tail_length]: default -1, print all reads of this contig. If this
 * 	value is specified, print only reads that overlaps either end of the
 * 	tail by [tail_length] bps
 */
void Contig::print (int contigID, const ReadIndexer& rIndexer,
	const Trimmer& myTrimmer, bool to_print_profile, bool to_print_reads,
	int tail_length) {
	/*
	 * print profile, consensus and read alignments
	 */
	//int batchSize = 1000000;
	ivec_t rIDs;
	get_rIDs (rIDs);
	std::vector<read_t> reads;
	rIndexer.get_trimmed_reads(std::back_inserter(reads), rIDs.begin(),
			rIDs.end(),	myTrimmer);

	std::set<ipair_t, comp_ipair> rIndices;
	for (int i = 0; i < (int) rIDs.size(); ++ i)
		rIndices.insert(ipair_t(rIDs[i], i));

	/* create profile and consensus */
	std::string consensus;
	ivec_t profile;
	generate_profile (consensus, profile, rIndices, reads);
	if (to_print_reads) print_reads (reads, tail_length);
	if(to_print_profile) print_profile (consensus, profile);

	/* consensus in fasta format */
	std::cout << ">" << contigID << "\t" << "length:" << consensus.length()
			<< "#reads:" << rIDs.size() << "\n";
	std::cout << consensus << "\n";
}

void Contig::print_profile(const std::string& consensus,
		const ivec_t& profile) const {
	for (int row = 0; row < 4; ++ row) {
		for (int col = 0; col < (int) profile.size()/4; ++ col) {
			std::cout << std::setw(5) << profile[4*col + row] << " ";
		}
		std::cout << "\n";
	}
	for (int i = 0; i < (int) consensus.length(); ++ i) {
		std::cout << std::setw(5) << consensus.at(i) << " ";
	}
	std::cout << "\n\n";
}

/* @brief	Reads in a contig is ordered by rID field. This function
 * 			instead sort reads by their distance to the beginning of the
 * 			contig, and the result is stored in a vector [dist_to_rID] of
 * 			<dist_to_beg, rID> pairs.
 */
void Contig::sort_by_dist (std::vector<ipair_t>& dist_to_rID) const {

	std::set<cread_t>::iterator it = db_.begin();
	for (; it != db_.end(); ++ it) {
		dist_to_rID.push_back(ipair_t(it->dist_to_beg, it->rID));
	}
	std::sort(dist_to_rID.begin(), dist_to_rID.end());
}

void Contig::print_reads (const std::vector<read_t>& reads,
		int tail_length) const {

	/* sort reads w.r.t. to dist_to_beg */
	std::vector<ipair_t> dist_to_rID;
	sort_by_dist (dist_to_rID);

	std::cout << "#reads: " << db_.size() << "\t";
	std::cout << "length: " << length_ << "\n\n";

	for (int i = 0; i < (int) dist_to_rID.size(); ++ i) {
		int rID = dist_to_rID[i].second;
		std::set<cread_t>::iterator it_r = db_.find(cread_t(rID));
		if (it_r == db_.end())
			abording("Contig:print_reads p1 sanity check failed");

		if (tail_length != -1) {
			int read_end = it_r->dist_to_beg + it_r->read_length +
						it_r->indel_length();
			if (it_r->dist_to_beg > tail_length ||
			   (read_end + tail_length >= length_)) continue;
		}

		for (int j = 0; j < it_r->dist_to_beg; ++ j) std::cout << "-";

		std::vector<read_t>::const_iterator
			iter = std::lower_bound(reads.begin(), reads.end(), read_t(rID));

		if (iter == reads.end())
			abording("Contig:print_reads p2 sanity check failed");

		std::string read = iter->read;
		std::string header = iter->header;

		if (!it_r->is_fwd_strand) read = xny::get_rvc_str (read);

		int ridx = 0;
		for (int j = 0; j < (int) it_r->indels.size(); ++ j) {
			int gap_start = it_r->indels[j].first;
			int gap_len = it_r->indels[j].second;
			std::cout << read.substr(ridx, gap_start - ridx);
			for (int g = 0; g < gap_len; ++ g) std::cout << "-";
			ridx = gap_start;
		}
		std::cout << read.substr(ridx, it_r->read_length - ridx) << "\t";

		if (it_r->is_fwd_strand) std::cout << "+\t";
		else std::cout << "-\t";
		std::cout << it_r->rID << "\t" << it_r->dist_to_beg
				<< "\t" << header << "\n";
	}
	std::cout << "\n";
}

/*
//[todo] remove [shingles] as the parameter, just for debugging purposes
void Contig::print_reads (const std::vector<shingle_t>& shingles) const {

	std::cout << "# reads: " << number_of_reads () << "\t";
	std::cout << "Length: " << length_ << "\n\n";

	for (int i = 0; i < number_of_reads(); ++ i) {
		for (int j = 0; j < creads_[i].dist_to_beg; ++ j) std::cout << "-";
		std::string read = creads_[i].seq;
		if (!creads_[i].is_fwd_strand) read = xny::get_rvc_str (read);

		std::cout << read.substr(0, creads_[i].gap_start);
		for (int gap = 0; gap < creads_[i].gap_length; ++ gap) std::cout << "-";
		std::cout << read.substr(creads_[i].gap_start,
			creads_[i].read_length - creads_[i].gap_start) << "\t";

		if (creads_[i].is_fwd_strand) std::cout << "+\t";
		else std::cout << "-\t";
		std::cout << creads_[i].dist_to_beg << "\t";

		int rID = creads_[i].rID;

		std::cout << rID << "\t" << shingles[rID].fwd_ss << "\t"
			<< (int) shingles[rID].fwd_min_pos << "\t" << shingles[rID].rv_ss << "\t"
			<< (int) shingles[rID].rv_min_pos << "\n";
	}
	std::cout << "\n";
}
*/

