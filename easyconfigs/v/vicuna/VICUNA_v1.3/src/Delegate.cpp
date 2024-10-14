//============================================================================
// Project     : Diversifier
// Name        : Delegate.cpp
// Author      : Xiao Yang
// Created on  : Oct 27, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#include "Delegate.h"

void Delegate::rvc () {
	int len = consensus.length();
	//std::swap(interval.first, interval.second);
	//interval.first = len - interval.first - 1;
	//interval.second = len - interval.second - 1;
	xny::rvc_str(consensus);
	for (int i = 0; i < len/2; ++ i) {
		/* swap a--c & g--t in ith and (len-i-1)th column, and
		 then these two columns, respectively */
		std::swap(profile[4*i + 0], profile[4*i + 3]);
		std::swap(profile[4*i + 1], profile[4*i + 2]);
		std::swap(profile[4*(len-i-1) + 0], profile[4*(len-i-1) + 3]);
		std::swap(profile[4*(len-i-1) + 1], profile[4*(len-i-1) + 2]);
		for (int j = 0; j < 4; ++ j)
			std::swap(profile[4*i + j], profile[4*(len-i-1) + j]);
	}
	if (len % 2 == 1) {
		int col = len/2;
		std::swap(profile[4*col + 0], profile[4*col + 3]);
		std::swap(profile[4*col + 1], profile[4*col + 2]);
	}

	if (lfv.size()) {
		std::reverse(lfv.begin(), lfv.end());
		for (int i = 0; i < lfv.size(); ++ i) {
			lfv[i].first = len - lfv[i].first - lfv[i].second;
		}
	}
} // rvc


/* @brief	Merge [rhs] to *this
 */
void Delegate::merge (const xny::align_t& aln, const Delegate& rhs) {
	std::string cons = "";
	int lhs_len = consensus.length(),
		rhs_len = rhs.length();
	ivec_t prfl;
	/* prefix */
	if (aln.s0_l <= aln.s1_l) { // lhs shorter
		cons += rhs.consensus.substr(0, aln.s1_l);
		prfl.insert(prfl.end(), rhs.profile.begin(),
					rhs.profile.begin() + 4*cons.length());
	} else {
		cons += consensus.substr(0, aln.s0_l);
		prfl.insert(prfl.end(), profile.begin(),
					profile.begin() + 4*cons.length());
	}

	/* alignment part */
	int col_0 = aln.s0_l,
		col_1 = aln.s1_l;
	for (int i = 0; i < (int) aln.seq.length(); ++ i) {
		ivec_t tmp(4, 0);
		int max_r, max_w;
		switch (aln.seq.at(i)) {
		case 'M':
			for (int r = 0; r < 4; ++ r)
				tmp[r] = profile[4*col_0+r] + rhs.profile[4*col_1+r];
			prfl.insert(prfl.end(), tmp.begin(), tmp.end());
			cons += consensus.at(col_0);
			++ col_0;
			++ col_1;
			break;
		case 'R':
			max_r = 0;
			max_w = -1;
			for (int r = 0; r < 4; ++ r) {
				tmp[r] = profile[4*col_0+r] + rhs.profile[4*col_1+r];
				if (tmp[r] > max_w) {
					max_r = r;
					max_w = tmp[r];
				}
			}
			prfl.insert(prfl.end(), tmp.begin(), tmp.end());
			cons += std::toupper(xny::bits2char(max_r));
			++ col_0;
			++ col_1;
			break;
		case 'I':
			prfl.insert(prfl.end(), rhs.profile.begin() + 4*col_1,
					rhs.profile.begin() + 4*col_1 + 4);
			cons += rhs.consensus.at(col_1);
			++ col_1;
			break;
		case 'D':
			prfl.insert(prfl.end(), profile.begin() + 4*col_0,
					profile.begin() + 4*col_0 + 4);
			cons += consensus.at(col_0);
			++ col_0;
			break;
		default:
			break;
		}
	}

	/* suffix */
	if (lhs_len - aln.s0_r <= rhs_len - aln.s1_r) { //lhs shorter
		if (aln.s1_r + 1 < rhs_len) {
			cons += rhs.consensus.substr(aln.s1_r+1,
					rhs_len - aln.s1_r - 1);
			prfl.insert(prfl.end(), rhs.profile.begin() +
					4 * (aln.s1_r + 1), rhs.profile.end());
		}

	} else {
		if (aln.s0_r + 1 < lhs_len) {
			cons += consensus.substr(aln.s0_r + 1,
				lhs_len - aln.s0_r - 1);
			prfl.insert(prfl.end(), profile.begin() +
				4 * (aln.s0_r + 1),	profile.end());
		}
	}

	consensus = cons;
	profile = prfl;
} // merge

/* @brief	Identify low frequent ([min_perc_pol]) variants
 * 			(could be sequencing errors) with constraint length
 * 			([max_variant_len]) for *this delegate. The results
 * 			are stored in private variable [lfv]
 * @method	First, identify all changing point (drop or increase of
 * 			coverage), Then, identify drop/increase pairs that satisfy
 * 			distance constraint. Nested pairs will be taken care of
 * 			this way.
 */
void Delegate::id_low_frq_var (int min_perc_pol, int max_variant_len) {
	// approach 1

	lfv.clear();
	int c_len = consensus.length();
	// position and type: drop or rise (left bracket or right bracket)
	ivec_t weights;
	std::vector<std::pair<int, char> > anchors;
	int pre_w = 0;
	for (int i = 0; i < 4; ++ i) pre_w += profile[i];
	pre_w = std::max(1, pre_w); // prevent dividing 0
	weights.push_back(pre_w);
	for (int col = 1; col < c_len; ++ col) {
		int cur_w = 0;
		for (int i = 0; i < 4; ++ i) cur_w += profile[4 * col + i];
		cur_w = std::max(1, cur_w);
		weights.push_back(cur_w);
		if (100 * cur_w / pre_w <= min_perc_pol) {
			// coverage dropping
			anchors.push_back(std::pair<int, char> (col, 'l'));
		} else if (100 * pre_w / cur_w <= min_perc_pol) {
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
				if (l > max_variant_len) break;
				if (anchors[i].second == 'l') {
					if (nested_r != 0) -- nested_r;
					else {
						int idx_e = anchors[index].first + 1;
						int idx_s = anchors[i].first - 1;
						int r_weight = weights[idx_e];
						int l_weight = weights[idx_s];
						bool low_frq = true;
						for (int j = idx_s + 2; j <= idx_e - 2; ++ j) {
							if (100 * weights[j]/l_weight > min_perc_pol
							|| 100 * weights[j]/r_weight > min_perc_pol){
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

	// approach 2
	/*
		lfv.clear();
		int c_len = consensus.length();
		int w_start = 1, w_one_bp_pre= INT_MAX, del_start = -1;
		for (int i = 0; i < 4; ++ i) w_start += profile[i];
		for (int col = 1; col < c_len; ++ col) {
			int w_end = 1;
			for (int i = 0; i < 4; ++ i) w_end += profile[4 * col + i];

			// check if cur is dropping
			if (100 * w_end/w_start < min_perc_pol) {
				if (del_start == -1) del_start = col;
				else if (col - del_start > max_variant_len) {
					w_start = w_end;
					del_start = -1;
				}
			} else {
				// check if pre stops dropping
				if (100 * w_one_bp_pre/ w_end < min_perc_pol && del_start != -1) {
					lfv.push_back(ipair_t (del_start, col - del_start));
				}
				w_start = w_end;
				del_start = -1;
			}
			w_one_bp_pre = w_end;
		}
	*/

	/*
	{ // debug printing
		if (lfv.size()) {
			std::string after_str;
			std::cout << "b4\n";
			for (int i = 0; i < lfv.size(); ++ i)
				std::cout << lfv[i].first <<"," << lfv[i].second << "\t";
			std::cout << "\n";
			print_profile ();
			std::cout << "\n>after\n";
			low_frq_var_free_consn (after_str);
			std::cout << after_str << "\n";
		}
	}
	*/
} // id_low_freq_variants

/* @brief	Get the consensus sequence [rc] after
 * 			removing low frq variants
 */
void Delegate::low_frq_var_free_consn (std::string& rcons) const {
	int c_len = consensus.length();
	rcons = "";
	if (lfv.size()) {
		int start = 0;
		for (int i = 0; i < lfv.size(); ++ i) {
			/*{ // debug
				std::cout << lfv[i].first << "," << lfv[i].second << "\t";
			}*/
			rcons += consensus.substr(start, lfv[i].first - start);
			start = lfv[i].first + lfv[i].second;
		}
		rcons += consensus.substr(start, c_len - start);
	} else {
		rcons = consensus;
	}
} // reduced_consensus

void Delegate::print_profile() {

	for (int i = 0; i < consensus.length(); ++ i) {
		std::cout << std::setw(6) << i << " ";
	}
	std::cout << "\n";

	for (int row = 0; row < 4; ++ row) {
		for (int col = 0; col < (int) profile.size()/4; ++ col) {
			std::cout << std::setw(6) << profile[4*col + row] << " ";
		}
		std::cout << "\n";
	}
	for (int i = 0; i < (int) consensus.length(); ++ i) {
		std::cout << std::setw(6) << consensus.at(i) << " ";
	}
	std::cout << "\n\n";
} // print_profile



/* @brief	Derive the reliable region of *this delegate,
 * 			where the start and the end profile columns satisfy
 * 			a) have the minimum weight of [min_prfl_col_weight]
 * 			b) the ratio b/t the weight of the dominant base to the total
 * 			weight is >= [min_cons_base_ratio]
 * 			c) start & end position is within [max_contig_overhang]
 * 			apart from the contig end.
 */
/*
void Delegate::get_reliable_region(int min_prfl_col_weight,
	int min_cons_base_ratio, int max_contig_overhang) {

	// compute the reliable region of the input contig by trimming
	 // off the beginning and the end using [min_consensus_base_ratio_]
	 // and   min_profile_col_weight_
	 //
	// forward direction
	int col_weight;
	int c_len = consensus.length();
	for (int col = 0; col < c_len; ++ col){
		col_weight = 0;
		for (int i = 0; i < 4; ++ i){
			col_weight += profile[4*col + i];
		}
		int consensus_row = xny::char2bits(consensus.at(col));
		if (consensus_row > 3)
			abording("Contiger::get_reliable_region: p0 SC failed");

		float cbase_weight = 100.0 * profile[4*col+consensus_row];
		if (col_weight >= min_prfl_col_weight &&
			cbase_weight/col_weight > min_cons_base_ratio) {
			ri.first = col;
			break;
		}
	}

	// reverse direction
	for (int col = c_len - 1; col >= 0; -- col){
		col_weight = 0;
		for (int i = 0; i < 4; ++ i)
			col_weight += profile[4*col + i];
		int consensus_row = xny::char2bits(consensus.at(col));
		if (consensus_row > 3)
			abording("Contiger::compute_individual_delegate: p1 SC failed");

		float cbase_weight = 100.0 * profile[4*col + consensus_row];
		if (col_weight >= min_prfl_col_weight &&
			cbase_weight/col_weight > min_cons_base_ratio) {
			ri.second = col;
			break;
		}
	}

	// max overhang criterion
	if (ri.first + 1 > max_contig_overhang) {
		ri.first = max_contig_overhang;
	}
	int remain_len = c_len - ri.second - 1;
	if (remain_len > max_contig_overhang) {
		ri.second = c_len - remain_len;
	}

	// reliable region is empty
	if (ri.second <= ri.first) {
		ri.second = ri.first = -1;
		return;
	}
} // get_reliable_region
*/
