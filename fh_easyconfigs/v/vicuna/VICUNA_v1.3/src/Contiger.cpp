//============================================================================
// Project     : Diversifier
// Name        : Contiger.cpp
// Author      : Xiao Yang
// Created on  : Aug 29, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#include "Contiger.h"
#include "time.h" // for random number generator
/*********************************************************************
 * @brief	The main function to build, merge & extend contigs
 *********************************************************************/
void Contiger::run (bool is_filter_on, const ivec_t& binned_rIDs,
		const Trimmer& myTrimmer, const ReadIndexer& rIndexer) {

	double timing = get_time();

	std::cout << "Building Contigs ...\n\n";

	std::cout << "\tcreate super shingles for " << binned_rIDs.size()
											   << " reads\n";

	std::vector<shingle_t> shingles;

	// [todo] is it necessary to clean [shingles] to remove entries
	// with readID = -1?
	super_shingling (shingles, binned_rIDs, rIndexer, myTrimmer);

    print_time("\t\t", timing);

	/* identify indices of [shingles] that share the same super-shingles */
	std::cout << "\tClustering via super-shingles\n";
	iivec_t common_shingle_indices;
	ss_clustering (common_shingle_indices, shingles);

    print_time("\t\t", timing);

	std::cout << "\tInit contigs\n";
	init_contig (common_shingle_indices, shingles, myTrimmer, rIndexer);

	print_time("\t\t", timing);


	//srand ( time(NULL) );
	//for (int i = 0; i < 3; ++ i) { // multiple fishing
	//	std::cout << "\t\tround: " << i << "\n";
		fishing (is_filter_on, myTrimmer, rIndexer);
	//}

	print_time("\t\t", timing);


	std::cout << "\tValidate contigs\n";
	validate_contig (myTrimmer, rIndexer);
	shingles.clear();
	common_shingle_indices.clear();

    print_time("\t\t", timing);

    std::cout << "\tMerge contigs via rID\n";
	/* union find structure to merge contigs */
	//merge_contig_via_rID_uf ();
	/* clustering approach to merge contigs */
	// merge_contig_via_rID ();
	merge_contig_via_rID_multi ();

    print_time("\t\t", timing);

     { // debug
    		int cnt = 0;
    		for (int i = 0 ; i < _contigs.size(); ++ i) {
    			if (_contigs[i].get_length() != 0) {
    				/*
    				std::cout << "\ncID:" << i <<"\t"
    						  << _contigs[i].get_length() << "\n";
    				std::set<cread_t>::iterator it = _contigs[i].get_db_begin();
    				for (;it != _contigs[i].get_db_end(); ++ it) {
    					std::cout << it->rID << "," << it->dist_to_beg << ","
    							<< (int) it->is_fwd_strand << ","
    							<< (int) it->read_length << "\n";
    				}
    				*/
    				++ cnt;
    			}
    		}
    		//std::cout << "cnt = " << cnt << "\n";
    	}

	std::cout << "\tValidate contigs\n";
	validate_contig (myTrimmer, rIndexer);

	/*
	{ // debug after merge_contig -- print out contigs
		reorganize_contigs ();
		std::cout << "\nTotal # Contigs: " << _contigs.size() << "\n";

		for (int i = 0; i < (int) _contigs.size(); ++ i) {
			iset_t contig_rIDs;
			_contigs[i].get_rIDs(contig_rIDs);
			if (contig_rIDs.count(303) || contig_rIDs.count(56616)) {
				std::cout << "contigID = " << i << "\trNum = " << contig_rIDs.size() << "\n";
				_contigs[i].print(i, rIndexer, myTrimmer, true, true);
			}
			//_contigs[i].print(i, rIndexer, myTrimmer, false, false);
		}
		exit(1);

		//std::cout << "\t\tTotal Reads involved in contigs: "
		  //          << contig_rIDs.size() << "\n";

		std::cout << "432\n\n";
		_contigs[432].print(432, rIndexer, myTrimmer, true, true);
		exit(1);
	}
	*/
	std::cout << "\tExtend contigs\n";
	extend_contig (myTrimmer, rIndexer);

    print_time("\t\t", timing);
} // run


/**********************************************************************
 * @brief	Generate super shingles for both fwd and rv strands of
 * 			each read if the read length is sufficient long
 **********************************************************************/
void Contiger::super_shingling (std::vector<shingle_t>& shingles,
		const ivec_t& rIDs, const ReadIndexer& rIndexer,
		const Trimmer& myTrimmer){

	// process reads in batches
	std::vector<read_t> reads;

	int numIDs = rIDs.size();
	int numBatches = numIDs/_batchsz;

	for (int i = 0; i < numBatches; ++ i) {
		rIndexer.getReads(std::back_inserter(reads),
				rIDs.begin() + i * _batchsz,
				rIDs.begin() + (i + 1) * _batchsz);

		batch_ss (shingles, reads, myTrimmer);

		reads.clear();
	}
	if (numIDs > numBatches * _batchsz) {
		rIndexer.getReads(std::back_inserter(reads),
				rIDs.begin() +  numBatches * _batchsz,	rIDs.end());

		batch_ss (shingles, reads, myTrimmer);
		reads.clear();
	}

} // super_shingling

void Contiger::batch_ss (std::vector<shingle_t>& shingles,
		const std::vector<read_t>& reads, const Trimmer& myTrimmer) {


	int size = reads.size();
	std::vector<shingle_t> localSS (size, shingle_t());
	std::string myRead;

	#pragma omp parallel for shared (localSS) private (myRead)

	for (int i = 0; i < size; ++ i) {
		myTrimmer.trim (myRead, reads[i]);
		localSS[i].rID = -1;
		if (generate_ss(localSS[i].fwd_ss, localSS[i].fwd_min_pos, myRead)
			&& generate_ss(localSS[i].rv_ss, localSS[i].rv_min_pos,
					xny::get_rvc_str(myRead))) {
			localSS[i].rID = reads[i].ID;
			localSS[i].is_paired = reads[i].paired;
		}
	}

	shingles.insert(shingles.end(), localSS.begin(), localSS.end());

} // batch_ss

bool Contiger::generate_ss (uint32_t& ss, uint8_t& common_shingle_pos,
							const std::string& read) {

	jaz::to_string ts;

	// stage 1
	int l = read.size();

	if (l - _w1 + 1 <= 0) return false;

	std::vector<std::pair<uint32_t, uint8_t> >
		sh1_pair(l - _w1 + 1, std::pair<uint32_t, uint8_t>());

	// jaz::murmur264 hash
	jaz::djb32 hash;

	for (int i = 0; i < l - _w1 + 1; ++i) {
		sh1_pair[i].first = hash(read.substr(i, _w1)); // std::string(read.begin() + i, read.begin() + i + _w1));
		sh1_pair[i].second = i;
	}

	std::sort(sh1_pair.begin(), sh1_pair.end(), comp_sh());

	std::vector<uint32_t> sh1(l - _w1 + 1, 0);
	for (int i = 0; i < l - _w1 + 1; ++ i) {
		sh1[i] = sh1_pair[i].first;
	}

	// stage 2
	l = sh1.size();

	if (l - _w2 + 1 <= 0) return false;

	std::vector<std::pair<uint32_t, uint8_t> >
		sh2_pair(l - _w2 + 1, std::pair<uint32_t, uint8_t>());

	//std::vector<uint32_t> sh2(l - _w2 + 1, 0);

    unsigned int v_sz = _w2 * sizeof(uint32_t);

    for (unsigned int i = 0; i < l - _w2 + 1; ++ i) {
    		const char* v = reinterpret_cast<const char*>(&sh1[i]);
    		//sh2[i] = hash(v, v_sz);
    		sh2_pair[i].first = hash(v, v_sz);
    		sh2_pair[i].second = i;
    }

    //ss = *std::min_element(sh2.begin(), sh2.end());
    std::sort(sh2_pair.begin(), sh2_pair.end(), comp_sh());
    ss = sh2_pair.front().first;
    int idx_min_ss = sh2_pair.front().second;
    common_shingle_pos = sh1_pair[idx_min_ss].second;
	return true;
} // generate_ss

void Contiger::ss_clustering (iivec_t& common_shingle_indices,
		const std::vector<shingle_t>& shingles) {

	std::map<uint32_t, ivec_t> ss_indices;
	std::map<uint32_t, ivec_t>::iterator it;
	for (int i = 0; i < (int) shingles.size(); ++ i) {
		if (shingles[i].rID == -1) continue;
		it = ss_indices.find(shingles[i].fwd_ss);
		if (it != ss_indices.end())	it->second.push_back(i);
		else ss_indices[shingles[i].fwd_ss] = ivec_t (1, i);

		it = ss_indices.find(shingles[i].rv_ss);
		if (it != ss_indices.end())	it->second.push_back(i);
		else ss_indices[shingles[i].rv_ss] = ivec_t (1, i);
	}

	// convert to vector structure
	for (it = ss_indices.begin(); it != ss_indices.end(); ++ it) {
		if (it->second.size() > 1) {
			common_shingle_indices.push_back(it->second);
		}
	}
} // ss_clustering

/************************************************************************
 * @brief	Initialize contigs based on min_hash
 ***********************************************************************/
void Contiger::init_contig(const iivec_t& common_shingle_indices,
		const std::vector<shingle_t>& shingles, const Trimmer& myTrimmer,
		const ReadIndexer& rIndexer){

	int counter = 0, startIdx = 0;
	iset_t batch_rIDs;
	std::vector<read_t> reads;

	iset_t all_rIDs; /* total no. reads involved in all contigs */

	for (int i = 0; i < (int) common_shingle_indices.size(); ++ i) {

		counter += common_shingle_indices[i].size();

		for (int j = 0; j < (int) common_shingle_indices[i].size(); ++ j) {
			batch_rIDs.insert(shingles[common_shingle_indices[i][j]].rID);
		}

		if (counter >= _batchsz) {

			//std::cout << "\t\tcounter = " << counter << "\n";

			rIndexer.get_trimmed_reads(std::back_inserter(reads),
					batch_rIDs.begin(), batch_rIDs.end(), myTrimmer);

			/*batch_build_contig (common_shingle_indices, startIdx, i,
					shingles, reads); */

			/* initialize current batch of contigs */
			for (int j = startIdx; j <= i; ++ j) {
				Contig candidate;
				candidate.init (common_shingle_indices[j], shingles, reads);
				_contigs.push_back(candidate);

				/* Stats for unique rIDs involved in contigs */
				iset_t tmp_rIDs;
				candidate.get_rIDs(tmp_rIDs);
				all_rIDs.insert(tmp_rIDs.begin(), tmp_rIDs.end());
			}

			startIdx = i + 1;

			batch_rIDs.clear();
			reads.clear();
			counter = 0;
		}
	} // for

	if (counter > 0) {

		//std::cout << "\t\tcounter = " << counter << "\n";

		rIndexer.get_trimmed_reads(std::back_inserter(reads),
				batch_rIDs.begin(), batch_rIDs.end(), myTrimmer);

		/* batch_build_contig (common_shingle_indices, startIdx,
				common_shingle_indices.size() - 1, shingles, reads); */

		/* initialize remaining contigs */
		int sz = common_shingle_indices.size();
		for (int j = startIdx; j <= sz - 1; ++ j) {
			Contig candidate;
			candidate.init (common_shingle_indices[j], shingles, reads);
			_contigs.push_back(candidate);
			/* Stats for unique rIDs involved in contigs */
			iset_t tmp_rIDs;
			candidate.get_rIDs(tmp_rIDs);
			all_rIDs.insert(tmp_rIDs.begin(), tmp_rIDs.end());
		}

	} // if (counter > 0)

	std::cout << "\t\tno. initial contigs: " << _contigs.size() << "\n"
			<< "\t\tno. reads involved: " <<  all_rIDs.size() << "\n";
} // init_contig

/***********************************************************************
/* @brief	Fishing out all remaining reads to merge into contigs or
 * 			form new contigs: if [filtered] = false, apply to raw input,
 * 			else, apply only to reads contained in contigs as well as
 * 			their paired reads.
 ***********************************************************************/
void Contiger::fishing (bool is_filter_on, const Trimmer& myTrimmer,
		const ReadIndexer& rIndexer) {

	std::vector<bool> seed_template;
	generate_seed (seed_template);
	
	{ // debug
		seed_template = std::vector<bool> (25, true);
		seed_template[2] = seed_template[5] = seed_template[11]
		    = seed_template[12] = seed_template[15] = seed_template[16]
		    = seed_template[17] = seed_template[18]= false;

	}
	
	
	/* print out seed */
	std::cout << "\t\tseed_template = ";
	for (int i = 0; i < (int) seed_template.size(); ++ i) {
		if (seed_template[i]) std::cout << "1";
		else std::cout << "0";
	}
	std::cout << "\n";

	/* seed -> (rID, rPos, dir) */
	std::map<uint64_t, spaced_seed_t> STable;

	/* generate rID_cID_map for all reads
	 * NOTE: although a read can belong to more than 1 existing
	 * contig, it doesn't matter which rID->cID get generated */
	imap_t rID_cID_map;
	generate_rID_cID_map (rID_cID_map, false);

	/* batch of input reads */
	ivec_t fish_rIDs;
	identify_reads_involved_in_fishing (fish_rIDs, rIndexer.size(),
			is_filter_on);

	std::vector<read_t> batch_reads;
	int num_batch = fish_rIDs.size() / _batchsz;

	for (int i = 0; i < num_batch; ++ i) {
		rIndexer.get_trimmed_reads(std::back_inserter(batch_reads),
				fish_rIDs.begin() + i * _batchsz,
				fish_rIDs.begin() + (i + 1) * _batchsz, myTrimmer);
		// process this batch
		batch_fishing (STable, rID_cID_map, seed_template, batch_reads);
		batch_reads.clear();
	}
	if (num_batch * _batchsz < fish_rIDs.size()) {
		rIndexer.get_trimmed_reads(std::back_inserter(batch_reads),
				fish_rIDs.begin() + num_batch * _batchsz,
				fish_rIDs.end(), myTrimmer);
		// process the remaining
		batch_fishing (STable, rID_cID_map, seed_template, batch_reads);
		batch_reads.clear();
	}
	fish_rIDs.clear();


	std::cout << "\t\tSTable size: " << STable.size() << "\n";

	/* remove singleton contigs */
	iset_t all_rIDs;
	int idx_empty_slot = -1;
	for (int i = 0; i < _contigs.size(); ++ i) {
		if (_contigs[i].number_of_reads() <= 1) {
			if (idx_empty_slot == -1) idx_empty_slot = i;
		} else {
			/* for stat purpose */
			iset_t tmp_rIDs;
			_contigs[i].get_rIDs(tmp_rIDs);
			all_rIDs.insert(tmp_rIDs.begin(), tmp_rIDs.end());
			_contigs[i].reset_db ();
			/* crunch empty contigs */
			if (idx_empty_slot != -1) {
				std::swap(_contigs[idx_empty_slot], _contigs[i]);
				++ idx_empty_slot;
			}
		}
	}
	
	// bug fixed by Xiao Yang Dec 11, 2012; when no contigs are singleton
	// idx_empty_slot remains to be -1, so no resizing shall be allowed
	if (idx_empty_slot != -1) {
		_contigs.resize(idx_empty_slot);
	}

	std::cout << "\t\tno. contigs: " << _contigs.size() << "\n"
			<< "\t\tno. reads involved: " <<  all_rIDs.size() << "\n";

} // fishing

/************************************************************************
 * @brief	Identify all rIDs involved in the fishing process,
 * 			if [filtered] = false, all input reads will participate,
 * 			otherwise, only reads already in some contig and their paired
 * 			ones will participate
 ***********************************************************************/
void Contiger::identify_reads_involved_in_fishing (ivec_t& rIDs,
		int input_read_cnt,	bool is_filter_on) {

	if (! is_filter_on) {
		rIDs.resize(input_read_cnt);
		for (int i = 0; i < input_read_cnt; ++ i) rIDs[i] = i;
	} else {
		std::set<cread_t>::const_iterator it_db;
		for (int i = 0; i < (int) _contigs.size(); ++ i) {
			it_db = _contigs[i].get_db_begin();
			for (; it_db != _contigs[i].get_db_end(); ++ it_db) {
				rIDs.push_back(it_db->rID);
				if (it_db->is_paired) {
					if (it_db->rID % 2 == 0) rIDs.push_back(it_db->rID + 1);
					else rIDs.push_back(it_db->rID - 1);
				}
			}
		} // for (int i
		std::sort (rIDs.begin(), rIDs.end());
		ivec_t::iterator it = std::unique_copy(rIDs.begin(), rIDs.end(),
				rIDs.begin());
		rIDs.resize(it - rIDs.begin());
	}
} // identify_reads_involved_in_fishing

/***********************************************************************
 * @brief	Generate seed for seeded-kmer matching
 * 			Window size = ( _w1 + _w2 ) * 1.5
 * 			no. '1's:  _w1 + _w2
 ***********************************************************************/
void Contiger::generate_seed (std::vector<bool>& seed) {

	int window_sz = (_w1 + _w2) * 1.5;
	seed = std::vector<bool> (window_sz, true);
	int zero_sz = window_sz - _w1 - _w2;
	ivec_t pos (window_sz - 2);
	for (int i = 0; i < (int) pos.size(); ++ i) pos[i] = i + 1;

	int times_shuffle = rand() % 37 + 1;
	for (int i = 0; i < times_shuffle; ++ i) {
		std::random_shuffle(pos.begin(), pos.end());
	}

	for (int i = 0; i < zero_sz; ++ i) seed[pos[i]] = false;

} // generate_seed


/* @brief
 *
 */
void Contiger::batch_fishing (std::map<uint64_t, spaced_seed_t>& STable,
		imap_t& rID_cID_map, const std::vector<bool>& seed_template,
		const std::vector<read_t>& reads) {

	int sz = reads.size();
	for (int i = 0; i < sz; ++ i) {

		/* query read in cread_t format */
		cread_t query_cread (reads[i].ID);
		query_cread.is_paired = reads[i].paired;
		query_cread.read_length = reads[i].read.length();

		/* [Si]: seeds of read i in the format of STable entries */
		std::map<uint64_t, spaced_seed_t> Si;

		/* get seeds of current read (both fwd & rvc strands) */
		std::vector<std::pair<uint64_t, int> > query_seeds;
		xny::get_bit_gapped_seeds (std::back_inserter(query_seeds),
				seed_template, reads[i].read, 2);

		/* search seeds in STable */
		bool found = false;
		for (int j = 0; j < (int) query_seeds.size(); ++ j) {

			uint64_t seed = query_seeds[j].first;

			spaced_seed_t Si_entry (reads[i].ID,
					(uint8_t) query_seeds[j].second, true);
			if (j % 2 != 0) { /* '-' strand */
				Si_entry.dir = false;
			}

			std::map<uint64_t, spaced_seed_t>::iterator
				it_hit_STable (STable.find(seed));

			if (it_hit_STable != STable.end()) { /* found */
				/* merge */
				//{ // debug	std::cout << "seed: " << seed << "\n";	}
				cluster_read_via_seed (query_cread, rID_cID_map,
						seed_template.size(), Si_entry, it_hit_STable);				found = true;
				break;
			} else Si[seed] = Si_entry;

		}

		if (!found) {

			/* add to database */
			STable.insert(Si.begin(), Si.end());

			/* create a new singleton contig && update rID_cID_map */

			Contig candidate;
			query_cread.dist_to_beg = 0;
			query_cread.is_fwd_strand = true;
			std::set<cread_t> db;
			db.insert(query_cread);
			candidate.set_db(db);
			candidate.set_length(query_cread.read_length);
			_contigs.push_back(candidate);

			rID_cID_map[query_cread.rID] = (int) _contigs.size() - 1;
		}

	}

} // batch_fishing

/**********************************************************************
 * @brief	Given query (R_q) and hit (R_h) reads that share a common
 * 			seed S', and contigs (C_q, C_h) they belong to, if R_q doesn't
 * 			belong to both C_q and C_h, then add R_q to C_h.
 *
 * @input	[query_cread]: query read in cread_t format
 * 			[rID_cID_map]: rID to cID mapping for all existing contigs
 * 			[seed_len]: length of seed template
 * 			[query_entry]: query read in spaced_seed_t format
 * 			[it_hit_STable]: iterator of hit in STable
 **********************************************************************/
void Contiger::cluster_read_via_seed (cread_t& query_cread,
	imap_t& rID_cID_map, int seed_len, const spaced_seed_t& query_entry,
	std::map<uint64_t, spaced_seed_t>::iterator it_hit_STable) {

	/* identify if query or hit rIDs exist in some contig */
	int q_rID = query_entry.rID,
		h_rID = it_hit_STable->second.rID;
	imap_t::iterator it_rc_q = rID_cID_map.find(q_rID);
	imap_t::iterator it_rc_h = rID_cID_map.find(h_rID);

	/* h_rID belongs to some contig */
	if (it_rc_h != rID_cID_map.end()) {

		int h_cID = it_rc_h->second;

		/* q_rID belongs to some contig and either query or hit contig
		 * contains both q_rID and h_rID */
		if (it_rc_q != rID_cID_map.end() &&
			((it_rc_q->second == it_rc_h->second) ||
			 (_contigs[h_cID].find_element(cread_t(q_rID)) !=
		 	  _contigs[h_cID].get_db_end())) ) return;

		//static int counter = 0;

		/* identify h_rID in the hit contig */
		std::set<cread_t>::iterator it_db =
			_contigs[h_cID].find_element(cread_t(h_rID));
		if (it_db == _contigs[h_cID].get_db_end()) {
			abording ("cluster_read_via_seeded_kmer: p1 SC failed");
		}

		/* adding [query_cread] to contig_[h_cID] */

		// 1) determine directionality
		query_cread.is_fwd_strand = it_db->is_fwd_strand;
		if (query_entry.dir != it_hit_STable->second.dir) {
			query_cread.is_fwd_strand =
					query_cread.is_fwd_strand ? false : true;
		}

		// 2) determine dist_to_beg
		int hit_offset = it_hit_STable->second.pos;
		int query_offset = query_entry.pos;

		if (it_db->is_fwd_strand) {  //h_rID.dir wrt contig +
			// seed on h_rID and q_rID: (+ -) or (- +)
			if (it_hit_STable->second.dir != query_entry.dir) {
				query_offset = (int) query_cread.read_length -
									query_offset - seed_len;
			}
		} else { // -
			hit_offset = (int) it_db->read_length - hit_offset - seed_len;
			// (+ +) or (- -)
			if (it_hit_STable->second.dir == query_entry.dir) {
				query_offset = (int) query_cread.read_length -
									query_offset - seed_len;
			}
		}

		/* NOTE: the dist_to_beg value for some reads can be negative,
		 * these values will be updated outside current function,
		 * as well as the length for each contig */
		query_cread.dist_to_beg = it_db->dist_to_beg +
										hit_offset - query_offset;

		// 3) apply addition
		_contigs[h_cID].insert(query_cread);

	} else {
		std::cout << "\tHit rID: " << h_rID << "belongs to no contig\n";
		abording ("in cluster_read_via_seed");
	}

} // cluster_read_via_seed

/**********************************************************************
 * @brief	Extend contigs using paired link structure;
 * 			iteratively choose the longest contig and extend with the
 *	 	    highest linked contig until no further extension can be done
 **********************************************************************/
void Contiger::extend_contig (const Trimmer& myTrimmer,
							  const ReadIndexer& rIndexer) {

	std::cout << "\t\tno. contigs to be extended: "
			  << _contigs.size() << "\n";
	/* sort contigs with decreasing length & clean */
	reorganize_contigs ();

	/* generate <rID, cID> map when rID is paired*/
	imap_t rID_cID_map;
	generate_rID_cID_map (rID_cID_map, true);

	/* [dgs]: delegates in the memory, which are updated if merging happens
	 * [add_dgs]: delegates to be added to [dgs]
	 * [ref_cID]: reference contig ID to be considered for merging */
	std::vector<Delegate> dgs, add_dgs;
	int ref_cID = 0;

	/* while more contigs to be processed */
	while (next_cID (ref_cID)) {


		{ // debug
	//		std::cout << "\n\t[[[extend contig "<< ref_cID << "]]]\n";

		}


		/* calculate neighbors of [ref_cID], to be stored in vec [ref_nbs],
		 * which is sorted in INCREASING # of links w/ [ref_cID] */
		std::vector<contig_nb_t> ref_nbs;
		get_contig_neighbors (ref_nbs, ref_cID, rID_cID_map);

		/*
		{ // debug print ref_nbs
			std::cout << ref_cID << ": nbs #=" << ref_nbs.size() << "\n";
			for (int i = 0; i < ref_nbs.size(); ++ i) {
				std::cout << ref_nbs[i].cID << "[" << ref_nbs[i].pp <<","
						<< ref_nbs[i].pm <<"," << ref_nbs[i].mp <<","
						<< ref_nbs[i].mm <<"]\n";
			}
			std::cout << "\n\n";
		}
		*/

		/* identify the index of [ref_cID] in [ref_nbs] */
		int ref_cID_nbs_idx = -1;
		if (! is_in_nbs (ref_cID_nbs_idx, ref_nbs, ref_cID)){
			ref_cID_nbs_idx = ref_nbs.size();
			ref_nbs.push_back(contig_nb_t(ref_cID));
		}

		/* identify delegates that are not yet in memory, add to [dgs],
		 * which is sorted by cID field */
		std::vector<contig_nb_t> add_nbs;
		compute_add_nbs (add_nbs, ref_nbs, dgs);
		if(add_nbs.size()){
			compute_delegates (add_dgs, add_nbs, myTrimmer, rIndexer);
			dgs.insert(dgs.end(), add_dgs.begin(), add_dgs.end());
			std::sort(dgs.begin(), dgs.end(), cmp_delegate());
			add_dgs.clear();
			add_nbs.clear();
		}

		/* remove ref_contig_id from nbs_ref */
		ref_nbs.erase(ref_nbs.begin() + ref_cID_nbs_idx);

		/*
		{	// debug print delegates
			for (int i = 0; i < dgs.size(); ++ i) {
				std::cout << ">" << dgs[i].cID << "\t len="
						<< dgs[i].consensus.length() << "\n"
						<< dgs[i].consensus << "\n\n";
			}
		}
		*/

		/* compare ref_contig_id with all of its neighbors */
		while (ref_nbs.size()) {

			int candidate_cID = ref_nbs.rbegin()->cID;

			/*
			{ // debug

				if (ref_cID == 64 && candidate_cID == 11) {
					std::cout << "here\n";
				}
			}
			*/

			if (ref_nbs.rbegin()->sum() < _min_contig_links) {
			//if (ref_nbs.rbegin()->sum() < 5) {
				ref_nbs.pop_back();
				continue;
			}
			/*
			{ // debug print
				std::cout << "nbs size = " << ref_nbs.size() << "\t";
				std::cout << "candidate_cID = " << candidate_cID << "\n";
			}
			*/

			bool rvc_cand = rvc_checker (*ref_nbs.rbegin());
			ref_nbs.pop_back();

			bool is_merged = false;
			/* candidate_cID can be merged to ref_cID ? */

			is_merged = merge_contig_via_overlap (ref_cID, candidate_cID,
												 rvc_cand, dgs);
			if (! is_merged) {
				/* rvc (candidate_cID) can be merged to ref_cID ? */
				rvc_cand = (rvc_cand == true) ? false: true;
				is_merged = merge_contig_via_overlap (ref_cID,
										candidate_cID, rvc_cand, dgs);
			}


			/*
			{ // debug print
				std::cout << "TRY MERGE " << ref_cID << "=" << candidate_cID <<"\t:";
				if (is_merged) std::cout << "succeed!\n\n";
				else std::cout << "failed!\n\n";
			}
			*/

			if (is_merged) {

				/*
				{ // debug
					if (ref_cID == 64 && candidate_cID == 60){
						std::cout << ">" << candidate_cID << ":\n";
						_contigs[ref_cID].print(ref_cID, rIndexer, myTrimmer, true, false);
					}
					if (ref_cID == 64 && candidate_cID == 11){
						std::cout << ">" << candidate_cID << ":\n";
						_contigs[ref_cID].print(ref_cID, rIndexer, myTrimmer, true, false);
						exit(1);
					}
				}
				*/

				/* identify nbs_cand & rm itself and ref_cID */
				std::vector<contig_nb_t> cand_nbs;
				get_contig_neighbors (cand_nbs, candidate_cID, rID_cID_map);

				int cand_cID_nbs_idx = -1;
				if (is_in_nbs (cand_cID_nbs_idx, cand_nbs, candidate_cID)){
					cand_nbs.erase(cand_nbs.begin() + cand_cID_nbs_idx);
				}
				if (is_in_nbs (ref_cID_nbs_idx, cand_nbs, ref_cID)) {
					cand_nbs.erase(cand_nbs.begin() + ref_cID_nbs_idx);
				}

				/* merge [cand_nbs] to [ref_nbs], which will be sorted
				 * in increasing order of total links */
				merge_nbs (ref_nbs, cand_nbs);

				/*
				{ // debug print ref_nbs
					std::cout << "after merging nbs, candidate contig: ";
					std::cout << candidate_cID << " nbs #= " << cand_nbs.size() << "\n";
					for (int i = 0; i < cand_nbs.size(); ++ i) {
						std::cout << cand_nbs[i].cID << "[" << cand_nbs[i].pp <<","
								<< cand_nbs[i].pm <<"," << cand_nbs[i].mp <<","
								<< cand_nbs[i].mm <<"]\n";
					}

					std::cout << "\n\n";

					std::cout << "overall ref_nbs: ";
					std::cout << ref_cID << ": nbs #= " << ref_nbs.size() << "\n";
					for (int i = 0; i < ref_nbs.size(); ++ i) {
						std::cout << ref_nbs[i].cID << "[" << ref_nbs[i].pp <<","
								<< ref_nbs[i].pm <<"," << ref_nbs[i].mp <<","
								<< ref_nbs[i].mm <<"]\n";
					}
					std::cout << "\n\n";

				}
				*/
				/* compute & include additional delegates if applicable */
				if (cand_nbs.size()) {

					compute_add_nbs (add_nbs, cand_nbs, dgs);

					if (add_nbs.size()){
						add_dgs = std::vector<Delegate> (add_nbs.size(), Delegate());
						compute_delegates (add_dgs, add_nbs, myTrimmer, rIndexer);
						dgs.insert(dgs.end(), add_dgs.begin(), add_dgs.end());
						std::sort(dgs.begin(), dgs.end(), cmp_delegate());
						add_nbs.clear();
						add_dgs.clear();
					}
				}
				/* update rcMap */
				update_rcMap (rID_cID_map, candidate_cID, ref_cID);

				/* clear candiate contig */
				_contigs[candidate_cID].clear();
			} // if (is_merged)

		} // while (ref_nbs.size())

		++ ref_cID;
	} //	while (next_cID () != -1)

	/* writing to output */
	output (dgs, myTrimmer, rIndexer);

} // extend_contig

/**********************************************************************
 * @brief	Output information of contigs, whose length is
 * 			>= [_min_output_contig_len]bp --
 * 			a) fasta format file containing contigs with
 * 			   lfv (low frequent variants)
 * 			b) corresponding alignment file of a)
 * 			c) fasta format file containing contigs w/o lfv
 *********************************************************************/
void Contiger::output (const std::vector<Delegate>& dgs,
		const Trimmer& myTrimmer, const ReadIndexer& rIndexer){

	std::string file_1 = _out_dir + "contig.lfv.fasta",
			file_2 = _out_dir + "contig.align",
			file_3 = _out_dir + "contig.fasta";
	std::ofstream o_1 (file_1.c_str()), o_2 (file_2.c_str()),
			o_3 (file_3.c_str());
	if (!o_1.good()) warning ("\tCan't create file: " + file_1);
	if (!o_2.good()) warning ("\tCan't create file: " + file_2);
	if (!o_3.good()) warning ("\tCan't create file: " + file_3);


	std::vector<ipair_t> len_dgIdx;
	for (int i = 0; i < (int) dgs.size(); ++ i) {
		if (dgs[i].length() >= _min_output_contig_len) {
			len_dgIdx.push_back(ipair_t(dgs[i].length(), i));
		}
	}
	std::sort(len_dgIdx.begin(), len_dgIdx.end(), comp_ipair());
	int sz = len_dgIdx.size();

	/* fasta format lfv inclusive contigs */
	if (o_1.good()) {
		for (int i = sz - 1; i >= 0; -- i) {
			int idx = len_dgIdx[i].second;
			o_1 << ">dg-" << sz - i - 1 << "\t"
				<< len_dgIdx[i].first << "\n"
				<< dgs[idx].get_consensus() << "\n";

			/*
			{ // debug

				if (len_dgIdx[i].first == 1792) {
					int cID = dgs[idx].get_cID();
					_contigs[cID].print(cID, rIndexer, myTrimmer, true, true, 150);
					exit(1);
				}
			}
			*/
		}
		o_1.close();
	}

	/* alignment information, lfv inclusive */
	if (o_2.good()) {
		int total_reads_involved = 0;
		for (int i = sz - 1; i >= 0; -- i) {
			int idx = len_dgIdx[i].second;
			int cID = dgs[idx].get_cID();
			_contigs[cID].print_alignment(o_2, sz - i - 1, rIndexer, myTrimmer);
			total_reads_involved += _contigs[cID].number_of_reads();
		}
		std::cout << "\t\tno. reads involved in final contigs: "
				  << total_reads_involved << "\n";
		o_2.close();
	}

	/* fasta format lfv free contigs */
	if (o_3.good()) {
		for (int i = sz - 1; i >= 0; -- i) {
			std::string s0;
			int idx = len_dgIdx[i].second;
			dgs[idx].low_frq_var_free_consn(s0);
			o_3 << ">dg-" << sz - i - 1 << "\t" << s0.length()
							<< "\tlfv_free\n" << s0 << "\n";
		}
		o_3.close();
	}

} // output

/* @brief	Given delegates [nbs], identify the ones that are not in
 * 			[dgs] and store these in [add_nbs]
 */
void Contiger::compute_add_nbs (std::vector<contig_nb_t>& add_nbs,
		const std::vector<contig_nb_t>& nbs,
		const std::vector<Delegate>& dgs) {

	if (dgs.size() == 0) add_nbs = nbs;
	else {
		for (int i = 0; i < (int) nbs.size(); ++ i) {
			if (std::binary_search(dgs.begin(), dgs.end(),
				Delegate(nbs[i].cID), cmp_delegate()) == false) {
				add_nbs.push_back(nbs[i]);
			}
		}
	}
} // compute_add_nbs


/**********************************************************************
 * @brief	Identify the next none zero length contig Idx in _contigs
 *			starting from [start_idx]
 * @para		[start_idx]: start position of the contig to be considered,
 * 			updated to be the first contig with none-zero length
 * @return	false if none found, otherwise, return the truex of the first
 * 			contig with none zero length
 **********************************************************************/
bool Contiger::next_cID(int& start_idx) {
	for (; start_idx < (int) _contigs.size(); ++ start_idx) {
		if (_contigs[start_idx].get_length() != 0) {
			return true;
		}
	}
	return false;
} // next_cID


/**********************************************************************
 * @brief	Repoint each rID in [src_cID] to [tgt_cID] in [rcMap]
 **********************************************************************/
void Contiger::update_rcMap (imap_t& rcMap, int src_cID, int tgt_cID){
	std::set<cread_t>::const_iterator
							it_db = _contigs[src_cID].get_db_begin();
	imap_t::iterator it_rc;
	for (; it_db != _contigs[src_cID].get_db_end(); ++ it_db) {
		it_rc = rcMap.find(it_db->rID);
		if (it_rc != rcMap.end()) {
			it_rc->second = tgt_cID;
		} else abording ("Contiger::update_rcMap SC failed");
	}
} // update_rcMap


/**********************************************************************
 * @brief	For each element of [rhs_nbs], IF it occurs in [lhs_nbs]
 * 			(identified by cID), accumulate each field of this element
 * 			to the one in [lhs_nbs] and rm it from [rhs_nbs], ELSE,
 * 			insert it to [lhs_nbs]. [lhs_nbs] is sorted in INCREASING
 * 			order of total links.
 * @param
 * 	[lhs_nbs]: target neighbor vector
 * 	[rhs_nbs]: source neighbor vector
 **********************************************************************/
void Contiger::merge_nbs (std::vector<contig_nb_t>& lhs_nbs,
		std::vector<contig_nb_t>& rhs_nbs){

	std::sort(lhs_nbs.begin(), lhs_nbs.end(), cmp_contig_nb_id());
	std::vector<contig_nb_t>::iterator it_lb;
	for (int i = 0; i < (int) rhs_nbs.size(); ++ i){
		it_lb = std::lower_bound(lhs_nbs.begin(), lhs_nbs.end(),
				rhs_nbs[i], cmp_contig_nb_id());
		if (it_lb != lhs_nbs.end()) {
			it_lb->mm += rhs_nbs[i].mm;
			it_lb->mp += rhs_nbs[i].mp;
			it_lb->pm += rhs_nbs[i].pm;
			it_lb->pp += rhs_nbs[i].pp;
			rhs_nbs.erase(rhs_nbs.begin() + i);
			-- i;
		}
	}
	lhs_nbs.insert(lhs_nbs.end(), rhs_nbs.begin(), rhs_nbs.end());
	std::sort(lhs_nbs.begin(), lhs_nbs.end(), cmp_contig_nb_lt());
} // merge_nbs


/**********************************************************************
 * @param
 * 	[nb]: info of paired links b/t *this to the reference contig
 * @return True if *this is likely to be in the rvc strand wrt the ref
 *********************************************************************/
bool Contiger::rvc_checker (const contig_nb_t& nb){
	int pp = nb.pp,	pm = nb.pm, mp = nb.mp, mm = nb.mm;
	if (((pp >= pm) && (pp >= mp)) || ((mm >= mp) && (mm >= pm))){
		return true;
	} else return false;
} // rvc_checker


/**********************************************************************
/* @param
 * 	[cID]: contig ID
 * 	[nbs]: paired link neighbors of some contig
 * @return	True if If [cID] is in [nbs], and its index [idx] in [nbs]
 **********************************************************************/
bool Contiger::is_in_nbs (int& idx, std::vector<contig_nb_t>& nbs,
		int cID){
	for (int i = (int) nbs.size() - 1; i >= 0; -- i) {
		if (nbs[i].cID == cID) {
			idx = i;
			return true;
		}
	}
	return false;
} // is_in_nbs


/**********************************************************************
 * @brief	Identify and remove ...ID... or ...DI... type of alignment
 * 			between s0 and s1, adjust aln coordinates accordingly.
 * @param
 * 	[aln]:	Alignment coordinates and sequence b/t [s0] and [s1]
 **********************************************************************/
void Contiger::close_alt_indels (xny::align_t& aln, const std::string& s0,
		const std::string& s1) {

	/* replace all ... ID ... structure into ...DI... */
	int pos = aln.seq.find("ID");
	while (pos != std::string::npos) {
		int D_cnt = 1, I_cnt = 1;
		for (int i = pos - 1; i >= 0; -- i) {
			if (aln.seq.at(i) == 'I') ++ I_cnt;
			else break;
		}
		for (int i = pos + 2; i < aln.seq.length(); ++ i) {
			if (aln.seq.at(i) == 'D') ++ D_cnt;
			else break;
		}
		if (D_cnt == I_cnt) {
			aln.seq.replace(pos - I_cnt + 1, 2 * D_cnt,
					std::string (D_cnt, 'M'));
		} else {
			aln.seq.replace(pos - I_cnt + 1, D_cnt + I_cnt,
					std::string (D_cnt, 'D') + std::string (I_cnt, 'I'));
		}
		pos = aln.seq.find("ID", pos + 1);
	}

	pos = aln.seq.find ("DI");
	while (pos != std::string::npos) {
		int D_cnt = 1, I_cnt = 1;
		int Ds = 0, Is = 0;
		for (int i = pos - 1; i >= 0; -- i) {
			if (aln.seq.at(i) == 'D') ++ D_cnt;
			else {
				/* calculating Ds and Is between beginning of alignment
				 * to current ..DI.. block */
				for (int j = i; j >= 0; -- j) {
					if (aln.seq.at(j) == 'I') ++ Is;
					else if (aln.seq.at(j) == 'D') ++ Ds;
				}
				break;
			}
		}
		for (int i = pos + 2; i < aln.seq.length(); ++ i) {
			if (aln.seq.at(i) == 'I') ++ I_cnt;
			else break;
		}

		int D_start = pos - D_cnt + 1;
		int I_start = pos + 1;
		if (D_cnt == I_cnt) {
			aln.seq.replace(D_start, 2 * D_cnt,	std::string (D_cnt, 'M'));
		} else {
			/* get the substring of s0 and s1 */
			std::string substr0 = s0.substr(aln.s0_l + D_start - Is, D_cnt),
					substr1 = s1.substr(aln.s1_l + I_start - Ds - D_cnt, I_cnt);
			if (D_cnt > I_cnt) {
				if (xny::hamming_dist(substr0.substr(0, I_cnt), substr1)
					< xny::hamming_dist(substr0.substr(D_cnt - I_cnt, I_cnt), substr1)) {
					aln.seq.replace(D_start, D_cnt + I_cnt,
						std::string (I_cnt, 'M') + std::string(D_cnt - I_cnt, 'D'));
				} else {
					aln.seq.replace(D_start, D_cnt + I_cnt,
						std::string(D_cnt - I_cnt, 'D') + std::string (I_cnt, 'M'));
				}
			} else {
				if (xny::hamming_dist(substr1.substr(0, D_cnt), substr0)
					< xny::hamming_dist(substr1.substr(I_cnt - D_cnt, D_cnt), substr0)) {
					aln.seq.replace(D_start, D_cnt + I_cnt,
						std::string (D_cnt, 'M') + std::string(I_cnt - D_cnt, 'I'));
				} else {
					aln.seq.replace(D_start, D_cnt + I_cnt,
						std::string(I_cnt - D_cnt, 'I') + std::string (D_cnt, 'M'));
				}
			}
		}

		pos = aln.seq.find("DI", pos + 1);
	}

} // close_alt_indels


/**********************************************************************
 * @brief	Try merge [cID1] to [cID0]
 * @para
 * 	[cID0]: target contig ID
 * 	[cID1]: source contig ID
 * 	[rvc1]: a suggestive flag that [cID1] may be in the rvc strand wrt
 * 			[cID0], however, if no valid overlap can be generated, the
 * 			other strand will be considered.
 * 	[delegates]: delegates for all contigs involved
 **********************************************************************/
bool Contiger::merge_contig_via_overlap (int cID0, int cID1, bool rvc1,
	std::vector<Delegate>& delegates) {

	bool debug = false;

	/* identify delegates for [cID1] and [cID0]*/
	std::vector<Delegate>::iterator
		it_dg0 = std::lower_bound(delegates.begin(), delegates.end(),
				Delegate(cID0), cmp_delegate()),
		it_dg1 = std::lower_bound(delegates.begin(), delegates.end(),
				Delegate(cID1), cmp_delegate());

	if (it_dg0 == delegates.end() || it_dg1 == delegates.end())
		abording ("Contiger::merge_contig_via_overlap: SC failed");

	// a
	//std::string cons_1 = it_dg1->get_consensus();
	//if (rvc1) xny::rvc_str(cons_1);

	/* remove low frq variants from consensus when applicable */
	///*
	//b
	std::string s0, s1;
	if (it_dg0->lfv_size()) it_dg0->low_frq_var_free_consn(s0);
	else s0 = it_dg0->get_consensus();
	if (it_dg1->lfv_size()) it_dg1->low_frq_var_free_consn(s1);
	else s1 = it_dg1->get_consensus();


	if (debug){// debug
		std::cout << "\n+++++Contigs:" << cID0 << " vs " << cID1 << "+++++\n";

		if (it_dg0->lfv_size()) {
			std::cout << ">b4 s0\t" << it_dg0->get_consensus().length()<<"\n";
			std::cout << it_dg0->get_consensus() << "\n";
			std::cout << ">after s0\t" << s0.length() << "\n" << s0 <<"\n\n";
		}
		if (it_dg1->lfv_size()) {
			std::cout << ">b4 s1\t" << it_dg1->get_consensus().length()<<"\n";
			std::cout << it_dg1->get_consensus() << "\n";
			std::cout << ">after s1\t" << s1.length() << "\n" << s1 <<"\n\n";
		}
	}

	if (rvc1) xny::rvc_str(s1);
	//*/


	/* kmer guided alignment */
	xny::align_t kframe;
	// a
	//xny::kmer_anchor_ps_aligner (kframe, it_dg0->get_consensus(), cons_1,
	//	seed_kmer_len_, min_contig_overlap_,	min_identity_, max_contig_overhang_);

	// b
	xny::kmer_anchor_ps_aligner (kframe, s0, s1, seed_kmer_len_,
		min_contig_overlap_,	min_identity_, max_contig_overhang_);
	/*** merging when applicable ***/
	if (!kframe.seq.length()) return false;
	else {

		/* Merge if kframe is close to one end of both contig */
		int overlap = std::min(kframe.s0_l, kframe.s1_l);
		int left_0 = kframe.s0_l, left_1 = kframe.s1_l;
		int right_0 = s0.length() - kframe.s0_r - 1,
			right_1 = s1.length() - kframe.s1_r - 1;
		overlap += std::min(right_0, right_1);
		overlap += std::max(kframe.s0_r - kframe.s0_l + 1,
						kframe.s1_r - kframe.s1_l + 1);

		int aln_len = std::min(kframe.s0_r - kframe.s0_l + 1,
				kframe.s1_r - kframe.s1_l + 1);

		if (100 * aln_len/overlap  >= 85 ||
			((kframe.s0_l <= max_contig_overhang_ ||
			  kframe.s1_l <= max_contig_overhang_) &&
			 (kframe.s0_r >= s0.length() - max_contig_overhang_ - 1||
			 kframe.s1_r >= s1.length() - max_contig_overhang_ - 1) &&
			 100 * aln_len/overlap  >= 75)
			 ) {
		/*
		if ( 100 * aln_len/overlap  >= 90 ||
			(left_0 >= left_1 && right_1 >= right_0 && left_1 <= 15) ||
			(left_1 >= left_0 && right_0 >= right_1 && right_1 <= 15) ||
			(left_0 >= left_1 && left_1 <= 15 && right_0 >= right_1 && right_1 <= 15) ||
			(left_0 <= left_1 && left_0 <= 15 && right_1 >= right_0 && right_0 <= 15)){
		*/
			/* reverse contig & delegate */
			if (rvc1) {
				_contigs[cID1].inverse_all_reads();
				it_dg1->rvc();
			}

			/* if [lfv] is not empty, convert s0 & s1 alignment
			 * to the original consensus alignment, adjust [kframe] */


			// b
			if (it_dg0->lfv_size()) { // s0


				if (debug){ // debug
					std::cout << "===S0===\n";
					for (std::vector<ipair_t>::const_iterator it = it_dg0->get_lfv_begin();
							it != it_dg0->get_lfv_end(); ++ it){
						std::cout << it->first << "," << it->second << "\t";
					}
					std::cout << "\n";
					std::cout << ">b4\t" << kframe.s0_l << "\t" << kframe.s0_r << "\n";

					if (cID0 == 124 && cID1 == 414) { // debug
						std::cout << kframe.seq << "\n\n";
					}
				}

				low_frq_var_rec (kframe.seq, kframe.s0_l, kframe.s0_r,
					it_dg0->get_lfv_begin(), it_dg0->get_lfv_end(), 'D');

				if (debug) {
					if (cID0 == 124 && cID1 == 414) { // debug
						std::cout << ">after\t" << kframe.s0_l << "\t" << kframe.s0_r << "\n";
						std::cout << kframe.seq << "\n\n";
					}
				}
			}

			if (it_dg1->lfv_size()) { //s1
				/*
				std::cout << "+++S1+++\n";
				for (std::vector<ipair_t>::const_iterator it = it_dg1->get_lfv_begin();
						it != it_dg1->get_lfv_end(); ++ it){
					std::cout << it->first << "," << it->second << "\t";
				}
				std::cout << "\n";

				std::cout << ">b4\t" << kframe.s1_l << "\t" << kframe.s1_r << "\n";
				std::cout << kframe.seq << "\n\n";
				*/

				low_frq_var_rec (kframe.seq, kframe.s1_l, kframe.s1_r,
					it_dg1->get_lfv_begin(), it_dg1->get_lfv_end(), 'I');
				/*
				std::cout << ">after\t" << kframe.s1_l << "\t" << kframe.s1_r << "\n";
				std::cout << kframe.seq << "\n\n";
				*/
			}

			/* close ID or DI gaps */
			// b
			if (it_dg0->lfv_size() || it_dg1->lfv_size()) {
				if (!rvc1) {
					close_alt_indels (kframe, it_dg0->get_consensus(),
							it_dg1->get_consensus());
				}
				else {
					close_alt_indels (kframe, it_dg0->get_consensus(),
							xny::get_rvc_str(it_dg1->get_consensus()));
				}

				/*
				std::cout << ">FURTHER\n";
				std::cout << kframe.seq << "\n\n";
				*/
			}

			if (debug) { // debug
				std::cout << "merging " << cID0 << ", " << cID1 << "\n";
			}

			_contigs[cID0].merge_via_overlap(kframe, _contigs[cID1]);

			/*
			{// debug
				std::cout << ">after length_ = \t"
						<< _contigs[cID0].get_length() << "\n";
			}
			*/
			/* e) merge [dg1] to [dg0] & remove [dg1] from delegates*/


			if (debug){// debug
				std::cout << "\n>b4\t" << it_dg0->length()
						<< "\n" <<  it_dg0->get_consensus() << "\n";
			}


			it_dg0->merge(kframe, *it_dg1);

			it_dg0->id_low_frq_var(min_perc_pol_, max_variant_len_);

			if (debug){// debug
				std::cout << ">after\t" << it_dg0->length()
						<< "\n" << it_dg0->get_consensus() << "\n\n"; // debug
			}

			delegates.erase(it_dg1);

			return true;

		} else {
			/* check if it is dove tail alignment or bad alignment */
			/*
			std::cout << "(" << s0.length() << "," << kframe.s0_l << "," << kframe.s0_r << ")\t";
			std::cout << "(" << s1.length() << "," << kframe.s1_l << "," << kframe.s1_r << ")\t";
			std::cout << "overhangs: \t" <<  std::min(kframe.s0_l, kframe.s1_l) << "," << std::min(s0.length() - kframe.s0_r - 1,
					s1.length() - kframe.s1_r - 1) << "\t";
			std::cout << "\t" << kframe.s0_r - kframe.s0_l + 1 << "\t" << overlap
						<< "\t" << 100 * (kframe.s0_r - kframe.s0_l + 1)/overlap << "\t";
			if (100 * (kframe.s0_r - kframe.s0_l + 1)/overlap  >= 85 || aln_len >= 60) std::cout << "********\n\n";
			else std::cout << "\n\n";
			*/

			return false; // ???
		}
	}

} // merge_contig_via_overlap


/* @brief	Let S' be the low-frq-variant-free sequence of contig S.
 * 			Given the [start] and the [end] positions of S' involved in
 * 			[aln], restore this alignment to the original S if
 * 			S' != S, i.e., there exist low-frq-variants (recorded as
 * 			<pos, len> pairs: [it_lfv_start], [it_lfv_end]) in S.
 * 			Character [ins] will be inserted to the [aln], and [aln_start]
 * 			[aln_end] will be updated accordingly
 */
void Contiger::low_frq_var_rec (std::string& aln, int& aln_start,
		int& aln_end, std::vector<ipair_t>::const_iterator it_lfv_start,
		std::vector<ipair_t>::const_iterator it_lfv_end, char ins) {

	int ins_cnt = 0; /* newly inserted Ds or Is */
	std::vector<ipair_t>::const_iterator	it = it_lfv_start;
	for (; it != it_lfv_end; ++ it) {
		if (it->first <= aln_start) {
			aln_start += it->second;
			aln_end += it->second;
		} else if (it->first <= aln_end){
			/* insert position wrt kframe.seq */
			int ins_pos = it->first - aln_start;
			for (int i = 0; i < aln.size(); ++ i) {
				if (aln.at(i) == 'M' || aln.at(i) == 'R' ||
						aln.at(i) == ins) {
					-- ins_pos;
					if (ins_pos < 0) {
						aln.insert(i, std::string(it->second, ins));
						break;
					}
				}
			} // for
			aln_end += it->second;
		} // if ... else if
	} // for

} // low_frq_var_rec



/*@brief		Generate delegates for all contigs by processing [_batchsz]
 * 			number of reads per batch
 */
void Contiger::compute_delegates (std::vector<Delegate>& delegates,
		const std::vector<contig_nb_t>& nbs, const Trimmer& myTrimmer,
		const ReadIndexer& rIndexer) {

	delegates.clear();
	delegates.resize(nbs.size());

	/* get all rIDs involved in the contigs in vector form so that
	 * the contig & reads ordering are preserved.
	 * [all_rIDs]: all rIDs involved in the contigs
	 * pref_sum[i]: num of reads up to the ith contig */
	ivec_t all_rIDs;
	ivec_t pref_sum (nbs.size(), 0);
	for (int i = 0; i < (int) nbs.size(); ++ i) {
		_contigs[nbs[i].cID].get_rIDs(all_rIDs);
		pref_sum[i] = all_rIDs.size();
	}

	/* load reads in batches */
	int total_rIDs = all_rIDs.size();
	int num_processed = 0;
	int nbs_idx_from, nbs_idx_to;

	const std::vector<Contig>* ptr_contigs = &_contigs;

	while (num_processed < total_rIDs) {

		iset_t batch_rIDs;
		int from = num_processed, to = num_processed + _batchsz;
		if (to < total_rIDs){
			batch_rIDs.insert(all_rIDs.begin() + from, all_rIDs.begin() + to);
			num_processed += _batchsz;
		} else {
			batch_rIDs.insert(all_rIDs.begin() + from, all_rIDs.end());
			to = (*pref_sum.rbegin());
			num_processed = total_rIDs;
		}

		/*identify the (start, end) of nbs chunks to be processed */
		ivec_t::iterator lb = std::lower_bound(pref_sum.begin(),
				pref_sum.end(), from);
		ivec_t::iterator ub = std::lower_bound(pref_sum.begin(),
				pref_sum.end(), to);

		if (lb == pref_sum.end() || ub == pref_sum.end()) {
			abording("Contiger::extend_contig p0 SC failed");
		} else {
			nbs_idx_from = (lb - pref_sum.begin());
			nbs_idx_to = (ub - pref_sum.begin());
		}
		std::vector<read_t> reads;
		rIndexer.get_trimmed_reads(std::back_inserter(reads),
				batch_rIDs.begin(), batch_rIDs.end(), myTrimmer);

		/*generate rIndices*/
		std::set<ipair_t, comp_ipair> rIndices;
		for (iset_t::iterator it = batch_rIDs.begin();
				it != batch_rIDs.end(); ++ it) {
			rIndices.insert(ipair_t(*it, rIndices.size()));
		}

		#pragma omp parallel for shared (ptr_contigs, rIndices, reads)

		for (int idx = nbs_idx_from; idx <= nbs_idx_to; ++ idx) {

			std::string consensus;
			ivec_t profile;

			/* option 1: weighted profile */
			(*ptr_contigs)[nbs[idx].cID].generate_weighted_profile(
					consensus, profile, rIndices, reads);

			/* option 2: normal profile */
			//(*ptr_contigs)[nbs[idx].cID].generate_profile(
			//		consensus, profile, rIndices, reads);

			delegates[idx].set_cID(nbs[idx].cID);
			delegates[idx].set_consensus(consensus);
			delegates[idx].set_profile (profile);
			delegates[idx].id_low_frq_var(min_perc_pol_,
											  max_variant_len_);
		}
	} // while (num_rIDs_processed < total_num_rIDs)

	std::sort(delegates.begin(), delegates.end(), cmp_delegate());
} // compute_delegates

/* @brief	For a given contig_[cID], generate paired link info b/t itself
 * 			and all linked contigs and singletons ([todo]),
 * @return	[nbs] sorted in the decreasing order of no. shared links
 */
void Contiger::get_contig_neighbors (std::vector<contig_nb_t>& nbs,
		int cID, const imap_t& rID_cID_map){

	// for debugging purposes
	//std::vector<ipair_t> pp, pm, mp, mm;

	std::vector<contig_nb_t>::iterator lower_bd;

	std::set<cread_t>::const_iterator it_db = _contigs[cID].get_db_begin(),
			it_db_nb;
	for (; it_db != _contigs[cID].get_db_end(); ++ it_db) {

		if (it_db->is_paired) {

			int neighbor_cID;
			imap_t::const_iterator it_rc;

			if (it_db->rID % 2 == 0) it_rc = rID_cID_map.find(it_db->rID + 1);
			else it_rc = rID_cID_map.find(it_db->rID - 1);

			/* the paired read is in another contig */
			if (it_rc != rID_cID_map.end()) {
				neighbor_cID = it_rc->second;

				it_db_nb = _contigs[neighbor_cID].find_element(
											cread_t(it_rc->first));

				if(it_db_nb == _contigs[neighbor_cID].get_db_end()) {

					std::cout << "rID: " << it_rc->first
							<< "\tnb cID: " << neighbor_cID << "\n";
					abording ("Contiger::get_contig_neighbors: "
							"p1 SC failed");
				}

				/* identify neighbor_cID in nbs, insert if not exist
				 * return existing or newly inserted iterator */
				lower_bd = std::lower_bound(nbs.begin(), nbs.end(),
						contig_nb_t(neighbor_cID), cmp_contig_nb_id());
				if (lower_bd != nbs.end()) {
					if (lower_bd->cID != neighbor_cID) { // not found
						lower_bd = nbs.insert(lower_bd,
								contig_nb_t(neighbor_cID));
					}
				} else {
					lower_bd = nbs.insert(nbs.end(),
							contig_nb_t(neighbor_cID));
				}

				if (it_db->is_fwd_strand) {
					if (it_db_nb->is_fwd_strand) {
						++ lower_bd->pp;
						//debug
						//if (neighbor_cID == cID)
						//pp.push_back(ipair_t(it_db->rID, it_rc->first));
					}
					else {
						++ lower_bd->pm;
						// debug
						//if (neighbor_cID == cID)
						//pm.push_back(ipair_t(it_db->rID, it_rc->first));
					}
				} else {
					if (it_db_nb->is_fwd_strand) {
						++ lower_bd->mp;
						// debug
						//if (neighbor_cID == cID)
						//mp.push_back(ipair_t(it_db->rID, it_rc->first));
					}
					else {
						++ lower_bd->mm;
						// debug
						//if (neighbor_cID == cID)
						//mm.push_back(ipair_t(it_db->rID, it_rc->first));
					}
				}
			}
			else { //[todo] the other pair doesn't belong to any contig

			}

		} // if (it_db->is_paired)

	} // for (; it_db != _contigs[cID].get_db_end(); ++ it_db)

	std::sort(nbs.begin(), nbs.end(), cmp_contig_nb_lt());

	// debug print out variety of links to itself
	/*
	std::cout << "\n------PP------\n";
	for (int i = 0; i < pp.size(); ++ i)
		std::cout << pp[i].first << "\t" << pp[i].second << "\n";
	std::cout << "\n------PM------\n";
	for (int i = 0; i < pm.size(); ++ i)
		std::cout << pm[i].first << "\t" << pm[i].second << "\n";
	std::cout << "\n------MP------\n";
	for (int i = 0; i < pp.size(); ++ i)
		std::cout << mp[i].first << "\t" << mp[i].second << "\n";
	std::cout << "\n------MM------\n";
	for (int i = 0; i < mm.size(); ++ i)
		std::cout << mm[i].first << "\t" << mm[i].second << "\n";
	*/
} // Contiger::get_contig_neighbors


/* @brief	Generate rID -> contigID map.
 * 			If [read_pair_only] = true, then only paired-reads
 * 			are considered, otherwise, all reads are considered.
 */
void Contiger::generate_rID_cID_map (imap_t& rID_cID_map,
		bool read_pair_only){

	std::set<cread_t>::const_iterator it_db;
	for (int i = 0; i < (int) _contigs.size(); ++ i) {
		it_db = _contigs[i].get_db_begin();
		for (; it_db != _contigs[i].get_db_end(); ++ it_db) {
			if (! read_pair_only)  rID_cID_map[it_db->rID] = i;
			else if (it_db->is_paired) {
				rID_cID_map[it_db->rID] = i;
			}
		}
	} // for (int i

} // Contiger::generate_rclink_map




/************************************************************************
 * @brief	Merge contig based on rID-> contig_indices and connected
 * 			component of clusters.
 * @method	Iterate:
 * 			1) Partition contigs into connected components.
 * 				a) select a reference contig C, identify all its
 * 				neighboring contigs	that share some common rID(s).
 * 				With C, they form a cluster.
 *		    		b) apply a) to the un-clustered contigs.
 *			2) For each cluster, use reference C and a base mark, update
 *				all other contigs (read offset, strand info), and merge
 *				with C.
 ***********************************************************************/
void Contiger::merge_contig_via_rID_multi () {

	int contig_num = _contigs.size();
	/* rID_list[i]: rIDs for _contigs[i] */
	std::vector<iset_t> rID_list (contig_num, iset_t());

	/* generate rID-> contig_idx_list map
	 * (a rID may belong to multiple contigs) */
	std::map<int, ivec_t> rID_cIDlist_map;
	std::map<int, ivec_t>::iterator it_rclist;
	//std::map<int, ipair_t> rID_cIDpair;
	//std::map<int, ipair_t>::iterator it_map;
	std::vector<ipair_t> rNum_cID; /* vector of <#reads, cID> pair */
	for (int idx = 0; idx < contig_num; ++ idx) {
		_contigs[idx].get_rIDs(rID_list[idx]);
		iset_t::iterator it = rID_list[idx].begin();
		for (; it != rID_list[idx].end(); ++ it){
			it_rclist = rID_cIDlist_map.find(*it);
			if (it_rclist != rID_cIDlist_map.end())
				it_rclist->second.push_back(idx);
			else rID_cIDlist_map[*it] = ivec_t(1, idx);
		}
		rNum_cID.push_back(ipair_t(rID_list[idx].size(), idx));
	}
	/* clean up entries where rID belongs to a unique contig */
	it_rclist = rID_cIDlist_map.begin();
	for (; it_rclist != rID_cIDlist_map.end(); ) {
		if (it_rclist->second.size() == 1){
			rID_cIDlist_map.erase(it_rclist ++);
		}
		else ++ it_rclist;
	}

	/* Handle contigs in an order of decreasing size (#reads) */
	std::sort(rNum_cID.begin(), rNum_cID.end(), comp_ipair());

	while (rID_cIDlist_map.size()) {

		/* partition of contigs */
		iivec_t cluster_of_cIDs;

		std::vector<bool> is_clustered (contig_num, false);
		for (int i = (int) rNum_cID.size() - 1; i >= 0; -- i) {

			/* [cluster]: stores all contigs that share common rIDs
			 * with [seed_cID] */
			ivec_t cluster;
			int seed_cID = rNum_cID[i].second;

			if (! is_clustered[seed_cID]) {

				cluster.push_back(seed_cID);
				is_clustered[seed_cID] = true;

				/* use every rID in the seed_cID to identify its neighbors */
				iset_t::iterator it = rID_list[seed_cID].begin();
				for (; it != rID_list[seed_cID].end(); ++ it) {
					it_rclist = rID_cIDlist_map.find(*it);
					if (it_rclist != rID_cIDlist_map.end()) {
						bool rc_flag = false;
						int rclist_sz = it_rclist->second.size();
						for (int j = 0; j < rclist_sz; ++ j) {
							int my_cID = it_rclist->second[j];
							if (my_cID != seed_cID) {
								if (! is_clustered[my_cID]){
									/* cluster any non-seed_cID contig */
									cluster.push_back(my_cID);
									is_clustered [my_cID] = true;
								}
 							} else rc_flag = true;
						}

						if (!rc_flag) {
							abording ("Contiger::merge_contig_via_rID "
									"SC failed");
						}

						rID_cIDlist_map.erase(it_rclist);

					} // if (it_rclist
				} // for (;

				if (cluster.size() > 1) {
					/* remove redundant cIDs from each [cluster] */
					std::sort (cluster.begin() + 1, cluster.end());
					ivec_t::iterator it = std::unique_copy(
							cluster.begin() + 1, cluster.end(),
							cluster.begin() + 1);
					cluster.resize(it - cluster.begin());
					cluster_of_cIDs.push_back(cluster);
				}

			} // if (! is_clustered[seed_cID])

		} // for (int i = (int) rNum_cID.size() - 1

		/* for omp purpose: redistribute the size of contig
		 * clusters to improve scalability during parallelization */
		/* std::random_shuffle(cluster_of_cIDs.begin(),
		  					   cluster_of_cIDs.end()); */

		merge_contigs_cc_multi (rID_cIDlist_map, rID_list, cluster_of_cIDs);

		/* update [rNum_cID] */
		rNum_cID.clear();
		for (int i = 0; i < (int) rID_list.size(); ++ i) {
			if (rID_list[i].size()) {
				rNum_cID.push_back(ipair_t (rID_list[i].size(), i));
			}
		}
		std::sort(rNum_cID.begin(), rNum_cID.end(), comp_ipair());
	} // while
} // merge_contig_via_rID_cc

/***********************************************************************
 * @brief	Merge every cluster of contigs
 * @method	For each [cluster], anchor all other contigs based on the ref
 * 			contig cluster[0] by updating constitutent read information.
 ***********************************************************************/
void Contiger::merge_contigs_cc_multi (std::map<int, ivec_t>& rID_cIDlist_map,
	std::vector<iset_t>& rID_list, const iivec_t& cluster_of_cIDs) {

	int sz = cluster_of_cIDs.size();

	/* this for loop can be parallelized by omp */
	for (int i = 0; i < sz; ++ i) {

		int seed_cID = cluster_of_cIDs[i][0];
		int seed_clen = _contigs[seed_cID].get_length();

		/* stores all read entries for the merged contig in vector form */
		std::vector<cread_t> all_entries;
		int min_dist = INT_MAX; // used for updating dist_to_beg

		/* handle contigs within each cluster */
		for (int j = 1; j < (int) cluster_of_cIDs[i].size(); ++ j) {

			std::vector<cread_t> add_entries;
			int common_rID = -1;

			int trg_cID = cluster_of_cIDs[i][j];
			int trg_clen = _contigs[trg_cID].get_length();

			/* check each entry in [trg_cID], identify 1 common rID with
			 * [seed_cID] and store all un-common rIDs to [add_entries],
			 * which will ultimated be merged to [seed_cID]
			 */

			std::set<cread_t>::iterator
						it_trg_db = _contigs[trg_cID].get_db_begin();

			int sd_dist_to_beg = -1, trg_dist_to_beg = -1;
			bool to_inverse = false;
			for (; it_trg_db != _contigs[trg_cID].get_db_end();
														++ it_trg_db) {

				iset_t::iterator it_sd_rID =
								rID_list[seed_cID].find(it_trg_db->rID);

				/* the rID is already included in rID_list[seed_cID] */
				if (it_sd_rID != rID_list[seed_cID].end()){

					/* the common_rID hasn't been identified yet */
					if (common_rID == -1){

						std::set<cread_t>::const_iterator
							it_sd_db = _contigs[seed_cID].find_element(
											cread_t(it_trg_db->rID));
						if (it_sd_db != _contigs[seed_cID].get_db_end()) {

							common_rID = it_trg_db->rID;
							to_inverse = (it_sd_db->is_fwd_strand ==
								it_trg_db->is_fwd_strand)? false : true;

							sd_dist_to_beg = it_sd_db->dist_to_beg;

							if (to_inverse) {
								trg_dist_to_beg = trg_clen -
										it_trg_db->dist_to_beg -
										it_trg_db->read_length;
							} else {
								trg_dist_to_beg = it_trg_db->dist_to_beg;
							}
						}
					}

				} else { /* new rID entry to be included */

					/* update rID_list for seed_cID */
					rID_list[seed_cID].insert(it_trg_db->rID);

					add_entries.push_back(*it_trg_db);
				}

				/* modify rID--> contigIdx mapping information */
				std::map<int, ivec_t>::iterator it_rclist =
						rID_cIDlist_map.find(it_trg_db->rID);

				if (it_rclist != rID_cIDlist_map.end()) {

					int rclist_sz = it_rclist->second.size();
					bool sc_flag = false;
					for (int j = 0; j < rclist_sz; ++ j) {
						if (it_rclist->second[j] == trg_cID) {
							it_rclist->second[j] = seed_cID;
							sc_flag = true;
						}
					}
					if (!sc_flag) {
						abording("Contiger::merge_contigs_cc SC p1 failed");
					}
				}

			} // for (it_trg_db

			if (common_rID == -1) {
				abording ("Contiger::merge_contig_cluster SC p2 failed: "
						"no common rID b/t seed_cID and trg_cID");
			}
			/* now update add_entries according to [to_inverse],
			 * [sd_dist_to_beg] [trg_dist_to_beg] */
			int offset = sd_dist_to_beg - trg_dist_to_beg;
			for (int l = 0; l < (int) add_entries.size(); ++ l) {
				if (to_inverse) {
					add_entries[l].is_fwd_strand =
							add_entries[l].is_fwd_strand? false:true;
					add_entries[l].dist_to_beg =
						trg_clen - add_entries[l].dist_to_beg -
						add_entries[l].read_length + offset;

				} else add_entries[l].dist_to_beg += offset;

				min_dist = std::min(min_dist, offset);
			} // for (l

			all_entries.insert(all_entries.end(), add_entries.begin(),
					add_entries.end());

			/* destroy [trg_cID] */
			_contigs[trg_cID].clear();
			rID_list[trg_cID].clear();
		} // for (j

		if (all_entries.size()) {
			/* update all_entries */
			if (min_dist < 0) {

				/* add the [seed_cID] to all_entries */
				all_entries.insert(all_entries.end(),
						_contigs[seed_cID].get_db_begin(),
						_contigs[seed_cID].get_db_end());

				_contigs[seed_cID].clear();

				for (int j = 0; j < (int) all_entries.size(); ++ j) {
					all_entries[j].dist_to_beg -= min_dist;
				}
			}

			_contigs[seed_cID].insert(all_entries.begin(), all_entries.end());
			_contigs[seed_cID].calculate_length();
		}

	} // for (i
} // merge_contigs_cc_multi



/* validation */
void Contiger::validate_contig (const Trimmer& myTrimmer,
							   const ReadIndexer& rIndexer) {

	std::vector<read_t> reads;
	iset_t batch_rIDs;
	imap_t rID_index;
	int startIdx = 0;
	int counter = 0;

	int size = _contigs.size();
	for (int i = 0; i < size; ++ i) {

		counter += _contigs[i].number_of_reads();
		_contigs[i].get_rIDs(batch_rIDs);
		if (counter >= _batchsz) {

			rIndexer.get_trimmed_reads(std::back_inserter(reads),
					batch_rIDs.begin(), batch_rIDs.end(), myTrimmer);

			int idx = 0;
			for (iset_t::iterator it = batch_rIDs.begin();
					it != batch_rIDs.end(); ++ it) {
				rID_index[*it] = idx;
				++ idx;
			}

			batch_validate_contig (startIdx, i, rID_index, reads,
									myTrimmer, rIndexer);

			startIdx = i + 1;
			reads.clear();
			batch_rIDs.clear();
			rID_index.clear();
			counter = 0;
		}
	}
	if (counter > 0) {

		rIndexer.get_trimmed_reads(std::back_inserter(reads),
				batch_rIDs.begin(), batch_rIDs.end(), myTrimmer);

		int idx = 0;
		for (iset_t::iterator it = batch_rIDs.begin();
				it != batch_rIDs.end(); ++ it) {
			rID_index[*it] = idx;
			++ idx;
		}

		batch_validate_contig (startIdx, size - 1, rID_index, reads,
				myTrimmer, rIndexer);
	}

	_contigs.erase(_contigs.begin(), _contigs.begin() + size);

	/* generate statistics */
	std::cout << "\t\tno. resulting contigs: " << _contigs.size() << "\n";
	iset_t all_rIDs;
	for (int i = 0; i < _contigs.size(); ++ i) {
		iset_t tmp_rIDs;
		_contigs[i].get_rIDs(tmp_rIDs);
		all_rIDs.insert(tmp_rIDs.begin(), tmp_rIDs.end());
	}
	std::cout << "\t\tno. reads involved: " << all_rIDs.size() << "\n";
}

/* @brief	Input [start] and [end] index of [contig_] to be considered
 * @param
 * 	[reads]: stores all possible read sequences for any rIDs involved in
 * 			the current batch.
 * 	[rID_index]: a map structure stores <rID, idx>, where idx is the index
 * 			of [reads] vector that stores the actual read sequence.
 */
void Contiger::batch_validate_contig (int start, int end,
		const imap_t& rID_index, const std::vector<read_t>& reads,
		const Trimmer& myTrimmer, const ReadIndexer& rIndexer) {

	std::vector<Contig> working_contigs (_contigs.begin() + start,
			_contigs.begin() + end + 1);
	//#pragma omp parallel shared (working_contigs) // a -- omp
	//{
		std::vector<Contig> localContigs;
		int size = working_contigs.size();
	//	#pragma omp for					  // a
		for (int i = 0; i < size; ++ i) {

			/*rIndices[i].first: rID; rIndices[i].second: index of [reads]*/
			std::set<ipair_t, comp_ipair> rIndices;

			std::set<cread_t>::const_iterator it_db = // a
					working_contigs[i].get_db_begin();  // a
			for (; it_db != working_contigs[i].get_db_end(); ++ it_db) { //a
				imap_t::const_iterator it_m = rID_index.find(it_db->rID);
				if (it_m != rID_index.end()) {
					rIndices.insert(ipair_t(it_m->first, it_m->second));
				} else abording ("Contiger::batch_validate_contig SC failed\n");
			}

			int init_num_reads = working_contigs[i].number_of_reads();

			std::set<ipair_t, comp_ipair> diverged;

			while (!working_contigs[i].validate (rIndices, reads, diverged,
					read_overhang_, percent_diverge_)) {

				localContigs.push_back(working_contigs[i].split(diverged));

				if (diverged.size() > 1) {
					rIndices = diverged;
					init_num_reads = working_contigs[i].number_of_reads();
					diverged.clear();
				}
				else break;
			}

			if ((int) diverged.size() >= init_num_reads - 1) continue;
			else { // candidate itself is a valid Contig if proper size
				if (working_contigs[i].number_of_reads() >= 2) {
					localContigs.push_back(working_contigs[i]);
				}
			}
		} // for (i =

		size = localContigs.size();
		for (int i = 0; i < size; ++ i)
			divide_contig_by_minoverlap (i, localContigs, rID_index, reads,
					myTrimmer, rIndexer);

	//	#pragma omp critical // a
	//	{
		_contigs.insert(_contigs.end(), localContigs.begin(), localContigs.end());
	//	}

	//} // #pragma omp parallel
}



/* @brief	Split the contig 'contigs[index]' into multiple whenever a
 * 			division partA and partB can be made such that their overlap
 * 			is < min_overlap (= _w1 + _w2, currently).
 * @para
 * 	[index]: index of contig in vector [contigs] to be considered
 * 	[contigs]: vector of contigs
 */
void Contiger::divide_contig_by_minoverlap (int index,
		std::vector<Contig>& contigs, const imap_t& rID_index,
		const std::vector<read_t>& reads,const Trimmer& myTrimmer,
		const ReadIndexer& rIndexer) {

	int min_overlap = _w1 + _w2;

	/* vector of (dist_to_beg, rID), which is sorted by dist_to_beg */
	std::vector<ipair_t> dist_to_rID;
	contigs[index].sort_by_dist(dist_to_rID);

	int contig_len = contigs[index].get_length();
	/* [pre_span]: stores the span of [tmp_contig] wrt contigs[index]
	 * [dist_to_beg]: the starting position of the [tmp_contig] wrt
	 * 				  the beginning of contigs[index].
	 * [tmp_contig]: stores the candidate contig
	 */
	int pre_span = 0, dist_to_beg = 0;

	Contig tmp_contig;
	ivec_t tmp_prfl = ivec_t (contig_len * 4, 0);
	std::string tmp_cons = std::string (contig_len, 'N');

	std::string read_seq;
	bool is_split = false;
	std::set<cread_t>::const_iterator it_db;

	for (int i = 0; i < (int) dist_to_rID.size(); ++ i) {

		int rID;

		rID = dist_to_rID[i].second;
		it_db = contigs[index].find_element(cread_t(rID));
		if (it_db == contigs[index].get_db_end()){
		  abording("Contiger::divide_contig_by_minoverlap p0 SC failed");
		}
		/* identify the actual read sequence */
		imap_t::const_iterator it_m = rID_index.find(rID);
		if (it_m != rID_index.end()) read_seq = reads[it_m->second].read;
		else abording ("Contiger::divide_contig_by_minoverlap SC p1 failed\n");

		/* initialize [tmp_contig] with the first read */
		if (i == 0) {
			pre_span = it_db->dist_to_beg + it_db->indel_length()
							+ it_db->read_length;
			dist_to_beg = it_db->dist_to_beg;

			tmp_contig.insert(*it_db);

			/* version 0: normal profile */
			tmp_contig.update_profile(tmp_cons, tmp_prfl, it_db, read_seq, 1);

			/* version 1: weighted profile */
			//tmp_contig.update_weighted_profile(tmp_cons, tmp_prfl,
			//								it_db, read_seq, 1);

		} else {

			/* the span of current read wrt contigs[index] */
			int cur_span = it_db->dist_to_beg + it_db->indel_length() +
						   it_db->read_length;

			int offset = pre_span - it_db->dist_to_beg;
			int overlap = std::min(offset, (int) read_seq.length());
			/* check the overlap b/t the current read w/ [tmp_contig] */

			std::string tmp_read = read_seq;
			if (!it_db->is_fwd_strand) xny::rvc_str(tmp_read);

			//if ((overlap >= (int) read_seq.length() * 0.6) ||
			//	(overlap >= min_overlap && (100.0 * xny::hamming_dist (
			//		tmp_cons.substr(	it_db->dist_to_beg, overlap),
			//		tmp_read.substr(0, overlap)) /overlap < 15))) { // good
			if (offset >= (int) read_seq.length() * 0.7 ||
				(offset >= min_contig_overlap_ &&
				  xny::hamming_dist(tmp_read.substr(0, overlap),
				    tmp_cons.substr(it_db->dist_to_beg, overlap))
					<= overlap * (100 - min_identity_) / 100)) { // good

				pre_span = std::max(pre_span, cur_span);
				tmp_contig.insert(*it_db);

				/* version 0: normal profile */
				tmp_contig.update_profile(tmp_cons, tmp_prfl, it_db, read_seq, 1);
				/* version 1: weighted profile */
				//tmp_contig.update_weighted_profile(tmp_cons, tmp_prfl,
				// 									it_db, read_seq, 1);
			} else { // split !


				/* calculate contig length */
				tmp_contig.set_length(pre_span - dist_to_beg);
				tmp_contig.update_dist_to_beg((-1) * dist_to_beg);

				/*
				{ // debug
					std::cout << "len = " <<  tmp_contig.get_length() << "\n";
					tmp_contig.print(0, rIndexer, myTrimmer, true, true);
				}
				*/
				contigs.push_back(tmp_contig);

				is_split = true;
				tmp_contig = Contig();
				pre_span = cur_span;
				dist_to_beg = it_db->dist_to_beg;
				tmp_contig.insert(*it_db);

				/* version 0: normal profile */
				tmp_contig.update_profile(tmp_cons, tmp_prfl, it_db, read_seq, 1);
				/* version 1: weighted profile */
				//tmp_contig.update_weighted_profile(tmp_cons, tmp_prfl,
				// 									it_db, read_seq, 1);

			}
		}
	}

	if (is_split) {
		tmp_contig.set_length(pre_span - dist_to_beg);
		tmp_contig.update_dist_to_beg((-1) * dist_to_beg);
		/*
		{ // debug
			std::cout << "len = " <<  tmp_contig.get_length() << "\n";
			tmp_contig.print(0, rIndexer, myTrimmer, true, true);
		}
		*/
		contigs[index] = tmp_contig;
	}
}


void Contiger::reorganize_contigs () {
	/* sort contig with increasing length */
	std::sort(_contigs.begin(), _contigs.end(), comp_contig_lt());

	/* identify the upper bound of the contigs with length 0 and resize */
	std::vector<Contig>::iterator up_bd = std::upper_bound(_contigs.begin(),
			_contigs.end(), Contig(), comp_contig_lt());

	_contigs.erase(_contigs.begin(), up_bd);

	/* sort contig with decreasing length */
	std::sort(_contigs.begin(), _contigs.end(), comp_contig_ge());
} // reorganize_contigs



/*************************************************************************
 * The following functions are replaced by new ones
/*************************************************************************/

/* split the contig into multiple if they are spaced
 *
 * store covered region of the contig in the format of
 * pairs of indices: (from, to) in the sorted order in [indices]
 *
 * initially [indices] is empty, for each read in the contig, update
 * the region it covered in [indices]
 */
void Contiger::divide_contig_by_n (int index, std::vector<Contig>& contigs){
	int contig_length = contigs[index].get_length();

	std::set<cread_t>::const_iterator it_db = contigs[index].get_db_begin();
	ivec_t indices;

	for (; it_db != contigs[index].get_db_end(); ++ it_db) {


		int from = it_db->dist_to_beg,
			to = it_db->dist_to_beg + it_db->indel_length() +
				 it_db->read_length - 1;

		ivec_t::iterator l_bd = std::lower_bound(indices.begin(),
												indices.end(), from);
		ivec_t::iterator u_bd = std::upper_bound(indices.begin(),
												indices.end(), to);

		if (l_bd != indices.end()) {

			int l_dist = (l_bd - indices.begin());

			if (u_bd != indices.end()) {
				int u_dist = (u_bd - indices.begin());

				indices.erase(l_bd, u_bd);

				if (l_dist % 2 == 0) {
					if (u_dist % 2 == 0) {
						indices.insert(indices.begin() + l_dist, to);
						indices.insert(indices.begin() + l_dist, from);
					} else {
						indices.insert(indices.begin() + l_dist, from);
					}
				} else {
					if (u_dist % 2 == 0) {
						indices.insert(indices.begin() + l_dist, to);
					}
				}
			} else {

				indices.erase(l_bd, u_bd);
				if (l_dist % 2 == 0) {
					indices.insert(indices.begin() + l_dist, to);
					indices.insert(indices.begin() + l_dist, from);
				} else indices.push_back(to);
			}

		} else {
			indices.push_back(from);
			indices.push_back(to);
		}

		if (indices.size() == 2 && indices[0] == 0 &&
				indices[1] == contig_length - 1) return;

	} // for (; it_db


	/* contig division */
	std::vector<Contig> split_contigs (indices.size()/2, Contig());
	it_db = contigs[index].get_db_begin();
	for (; it_db != contigs[index].get_db_end(); ++ it_db) {
		int from = it_db->dist_to_beg,
			to = it_db->dist_to_beg+it_db->indel_length()+it_db->read_length-1;
		ivec_t::iterator l_bd = std::lower_bound(indices.begin(), indices.end(), from);
		int l_dist = (l_bd - indices.begin());
		if (l_dist % 2 == 0) {
			if (from >= indices[l_dist] && to <= indices[l_dist + 1]) {
				split_contigs[l_dist/2].insert(*it_db);
			} else abording("Contiger::divide_contig_by_n: p1 Sanity check failed");

		} else {
			if (from >= indices[l_dist - 1] && to <= indices[l_dist]) {
				split_contigs[l_dist/2].insert(*it_db);
			} else abording("Contiger::divide_contig_by_n: p2 Sanity check failed");

		}
	}
	for (int i = 0; i < (int) split_contigs.size(); ++ i) {
		split_contigs[i].update_dist_to_beg(-1 * indices[i*2]);
		split_contigs[i].set_length(indices[i*2+1] - indices[i*2] + 1);
	}

	contigs[index] = split_contigs[0];
	contigs.insert(contigs.end(), split_contigs.begin() + 1, split_contigs.end());

} // Contiger::divide_contig_by_n

/* @brief	Split the contig 'contigs[index]' into multiple whenever a
 * 			division partA and partB can be made such that their overlap
 * 			is < min_overlap (= _w1 + _w2, currently).
 * @para
 * 	[index]: index of contig in vector [contigs] to be considered
 * 	[contigs]: vector of contigs
 */
void Contiger::divide_contig_by_minoverlap_0 (int index,
		std::vector<Contig>& contigs, const imap_t& rID_index,
		const std::vector<read_t>& reads) {
	int min_overlap = _w1 + _w2;

	/* vector of (dist_to_beg, rID), which is sorted by dist_to_beg */
	std::vector<ipair_t> dist_to_rID;
	contigs[index].sort_by_dist(dist_to_rID);

	int contig_len = contigs[index].get_length();
	/* [pre_span]: stores the span of [tmp_contig] wrt contigs[index]
	 * [dist_to_beg]: the starting position of the [tmp_contig] wrt
	 * 				  the beginning of contigs[index].
	 * [tmp_contig]: stores the candidate contig
	 */
	int pre_span = 0, dist_to_beg = 0;
	Contig tmp_contig;
	ivec_t tmp_prfl = ivec_t (contig_len * 4, 0);
	std::string tmp_cons = std::string (contig_len, 'N');
	std::string read_seq;
	bool is_split = false;

	std::set<cread_t>::const_iterator it_db;

	for (int i = 0; i < (int) dist_to_rID.size(); ++ i) {

		int rID;

		/* initialize [tmp_contig] with the first read */
		if (i == 0) {
			rID = dist_to_rID[i].second;
			it_db = contigs[index].find_element(cread_t(rID));
			if (it_db == contigs[index].get_db_end()){
			  abording("Contiger::divide_contig_by_minoverlap p1 SC failed");
			}
			pre_span = it_db->dist_to_beg + it_db->indel_length()
							+ it_db->read_length;
			dist_to_beg = it_db->dist_to_beg;

			tmp_contig.insert(*it_db);

			/* identify the actual read sequence */
			imap_t::const_iterator it_m = rID_index.find(rID);
			if (it_m != rID_index.end()) {
				read_seq = reads[it_m->second].read;
				tmp_contig.update_weighted_profile(tmp_cons, tmp_prfl,
												it_db, read_seq, 1);
			} else {
				abording ("Contiger::divide_contig_by_minoverlap SC p0 failed\n");
			}

			continue;
		}

		rID = dist_to_rID[i].second;
		it_db = contigs[index].find_element(cread_t(rID));
		if (it_db == contigs[index].get_db_end()) {
			std::cout << "rID = " << rID << "\n";
			for (int x = 0; x < dist_to_rID.size(); ++ x) {
				std::cout << dist_to_rID[x].first << " -- "
						<< dist_to_rID[x].second << "\n";
			}
			abording("Contiger::divide_contig_by_minoverlap p2 SC failed");
		}

		/* the span of current read wrt contigs[index] */
		int cur_span = it_db->dist_to_beg + it_db->indel_length() +
					   it_db->read_length;

		/* identify the actual read sequence */
		imap_t::const_iterator it_m = rID_index.find(rID);
		if (it_m != rID_index.end()) read_seq = reads[it_m->second].read;
		else abording ("Contiger::divide_contig_by_minoverlap SC p1 failed\n");

		int overlap = pre_span - it_db->dist_to_beg;
		/* check the overlap b/t the current read w/ [tmp_contig] */
		if (overlap >= min_overlap &&
			(100.0 * xny::hamming_dist (tmp_cons.substr(
				tmp_cons.length() - overlap, overlap),
				read_seq.substr(0, overlap)) /overlap < 10)) { // good
			pre_span = std::max(pre_span, cur_span);
			tmp_contig.insert(*it_db);
			tmp_contig.update_weighted_profile(tmp_cons, tmp_prfl,
											it_db, read_seq, 1);
		} else { // split !
			/* calculate contig length */
			tmp_contig.set_length(pre_span - dist_to_beg);
			tmp_contig.update_dist_to_beg((-1) * dist_to_beg);
			contigs.push_back(tmp_contig);

			is_split = true;
			tmp_contig = Contig();
			pre_span = cur_span;
			dist_to_beg = it_db->dist_to_beg;
			tmp_contig.insert(*it_db);

			ivec_t tmp_prfl = ivec_t (contig_len * 4, 0);
			std::string tmp_cons = std::string (contig_len, 'N');
			tmp_contig.update_weighted_profile(tmp_cons, tmp_prfl,
												it_db, read_seq, 1);

		}
	}

	if (is_split) {
		tmp_contig.set_length(pre_span - dist_to_beg);
		tmp_contig.update_dist_to_beg((-1) * dist_to_beg);
		contigs[index] = tmp_contig;
	}
}


/* @brief 	Merge contigs based on rID-> contigs_indices, where
 * 			a union-find method is used (this method is replace by
 * 			merge_contig_via_rID_cc() to speed up the process.
 */
void Contiger::merge_contig_via_rID_uf () {

	/* union find structure to guide merging of contigs
	 * [todo] This can be parallelized by partitioning the task */

	/* initialize to the number of contigs */
	ivec_t ufvec (_contigs.size(), 0);
	for (int i = 0; i < (int) ufvec.size(); ++ i) ufvec[i] = i;

	/* generate rID-> contig_indices (ipair_t), since each
	   rID belongs to at most 2 contigs */
	std::map<int, ipair_t> rID2contigIdx;
	std::map<int, ipair_t>::iterator it_m;
	for (int idx = 0; idx < (int) _contigs.size(); ++ idx) {
		iset_t rIDs;
		_contigs[idx].get_rIDs(rIDs);
		for (iset_t::iterator it_s = rIDs.begin(); it_s != rIDs.end(); ++ it_s){
			it_m = rID2contigIdx.find(*it_s);
			if (it_m != rID2contigIdx.end()) it_m->second.second = idx;
			else rID2contigIdx[*it_s] = ipair_t(idx, -1);
		}
	}

	//int debug_counter = 0;
	for (it_m = rID2contigIdx.begin(); it_m != rID2contigIdx.end(); ++ it_m) {

		/*
		if (debug_counter % 1000 == 0)
			std::cout << "counter = " << debug_counter << "\n";
		debug_counter ++ ;
		*/
		int common_rID = it_m->first;
		if (it_m->second.second != -1) {
			int contig_idx1 = it_m->second.first;
			int root_1 = find (contig_idx1, ufvec);
			int contig_idx2 = it_m->second.second;
			int root_2 = find (contig_idx2, ufvec);
			if (root_1 != root_2) {
				/* merge the contig with less # reads to the one with larger */
				if (_contigs[root_1].number_of_reads() >=
						_contigs[root_2].number_of_reads()) {
					if (_contigs[root_1].merge_via_rID(_contigs[root_2],
							common_rID)) {
						ufvec[root_2] = root_1;
					}
				} else if (_contigs[root_2].merge_via_rID(_contigs[root_1],
						common_rID)) {
					ufvec[root_1] = root_2;
				}
			} // if (root_1 != root_2)
		} // if
	} // for it_m

} // Contiger::merge_contig


/* @brief	Merge contig based on rID-> contig_indices and connected
 * 			component of clusters.
 * @method	Iterate:
 * 			1) Partition contigs into connected components.
 * 				a) select a reference contig C, identify all its
 * 				neighboring contigs	that share some common rID(s).
 * 				With C, they form a cluster.
 *		    		b) apply a) to the un-clustered contigs.
 *			2) For each cluster, use reference C and a base mark, update
 *				all other contigs (read offset, strand info), and merge
 *				with C.
 */
void Contiger::merge_contig_via_rID () {

	int contig_num = _contigs.size();
	/* rID_list[i]: rIDs for _contigs[i] */
	std::vector<iset_t> rID_list (contig_num, iset_t());

	/* generate rID-> contig_idx_pair (a rID belongs to max 2 contigs) */
	std::map<int, ipair_t> rID_cIDpair;
	std::map<int, ipair_t>::iterator it_map;
	std::vector<ipair_t> size_cID; /* pairs of (#reads, cID) */
	for (int idx = 0; idx < contig_num; ++ idx) {
		_contigs[idx].get_rIDs(rID_list[idx]);
		iset_t::iterator it = rID_list[idx].begin();
		for (; it != rID_list[idx].end(); ++ it){
			it_map = rID_cIDpair.find(*it);
			if (it_map != rID_cIDpair.end()) it_map->second.second = idx;
			else rID_cIDpair[*it] = ipair_t(idx, -1);
		}
		size_cID.push_back(ipair_t(rID_list[idx].size(), idx));
	}
	/* clean up rIDs that belong to a unique contig */
	for (it_map = rID_cIDpair.begin(); it_map != rID_cIDpair.end(); ) {
		if (it_map->second.second == -1) rID_cIDpair.erase(it_map ++);
		else ++ it_map;
	}

	/* Handle contigs in an order of decreasing size (#reads) */
	std::sort(size_cID.begin(), size_cID.end(), comp_ipair());

	while (rID_cIDpair.size()) {

		/* partition of contigs */
		iivec_t cluster_of_cIDs;

		std::vector<bool> is_clustered (contig_num, false);
		for (int i = (int) size_cID.size() - 1; i >= 0; -- i) {

			/* cluster of contig IDs, cluster[0] is the ref contig for
			 * initiating the clustering */
			ivec_t cluster;
			int seed_cID = size_cID[i].second;

			if (! is_clustered[seed_cID]) {

				cluster.push_back(seed_cID);
				is_clustered[seed_cID] = true;

				/* use every rID in the seed_cID to identify its neighbors */
				iset_t::iterator it = rID_list[seed_cID].begin();
				for (; it != rID_list[seed_cID].end(); ++ it) {
					it_map = rID_cIDpair.find(*it);
					if (it_map != rID_cIDpair.end()) {

						if (it_map->second.first == seed_cID) {
							if (it_map->second.second == seed_cID) {
								rID_cIDpair.erase(it_map);
							} else if (! is_clustered[it_map->second.second]) {
								cluster.push_back(it_map->second.second);
								is_clustered[it_map->second.second] = true;
								rID_cIDpair.erase(it_map);
							}
						} else if (it_map->second.second == seed_cID) {
							if (! is_clustered[it_map->second.first]) {
								cluster.push_back(it_map->second.first);
								is_clustered[it_map->second.first] = true;
								rID_cIDpair.erase(it_map);
							}
						} else {
							abording ("Contiger::merge_contig_via_rID SC failed");
						}
					} // if (it_map
				} // for (;

				if (cluster.size() > 1) {
					/* remove redundant cIDs from each [cluster] */
					std::sort (cluster.begin() + 1, cluster.end());
					ivec_t::iterator it = std::unique_copy(
							cluster.begin() + 1, cluster.end(),
							cluster.begin() + 1);
					cluster.resize(it - cluster.begin());
					cluster_of_cIDs.push_back(cluster);
				}
			} // if (! is_clustered[seed_cID])
		} // for (int i = (int) size_cID.size() - 1; i >= 0; -- i)

		/* redistribute the size of contig clusters (improve parallelization) */
		std::random_shuffle(cluster_of_cIDs.begin(), cluster_of_cIDs.end());

		merge_contigs_cc (rID_cIDpair, rID_list, cluster_of_cIDs);

		/* update [size_cID] */
		size_cID.clear();
		for (int i = 0; i < (int) rID_list.size(); ++ i) {
			if (rID_list[i].size()) {
				size_cID.push_back(ipair_t (rID_list[i].size(), i));
			}
		}
		std::sort(size_cID.begin(), size_cID.end(), comp_ipair());
	} // while
} // merge_contig_via_rID

/* @brief	Merge every cluster of contigs
 * @method	For each [cluster], anchor all other contigs based on the ref
 * 			contig cluster[0] by updating constitutent read information.
 */
void Contiger::merge_contigs_cc (std::map<int, ipair_t>& rID_cIDpair,
	std::vector<iset_t>& rID_list, const iivec_t& cluster_of_cIDs) {

	int sz = cluster_of_cIDs.size();

	/* this for loop can be parallelized by omp */
	for (int i = 0; i < sz; ++ i) {

		int seed_cID = cluster_of_cIDs[i][0];
		int seed_clen = _contigs[seed_cID].get_length();

		/* stores all read entries for the merged contig in vector form */
		std::vector<cread_t> all_entries;
		int min_dist = INT_MAX;

		/* handle contigs within each cluster */
		for (int j = 1; j < (int) cluster_of_cIDs[i].size(); ++ j) {

			std::vector<cread_t> add_entries;
			int common_rID = -1;

			int trg_cID = cluster_of_cIDs[i][j];
			int trg_clen = _contigs[trg_cID].get_length();

			/* check each entry in [trg_cID], identify 1 common rID with
			 * [seed_cID] and store all un-common rIDs to [add_entries],
			 * which will ultimated be merged to [seed_cID]
			 */

			std::set<cread_t>::iterator
						it_trg_db = _contigs[trg_cID].get_db_begin();

			int sd_dist_to_beg = -1, trg_dist_to_beg = -1;
			bool to_inverse = false;
			for (; it_trg_db != _contigs[trg_cID].get_db_end();
														++ it_trg_db) {

				iset_t::iterator it_sd_rID =
								rID_list[seed_cID].find(it_trg_db->rID);

				/* the rID is already included in rID_list[seed_cID] */
				if (it_sd_rID != rID_list[seed_cID].end()){

					/* the common_rID hasn't been identified yet */
					if (common_rID == -1){

						std::set<cread_t>::const_iterator
							it_sd_db = _contigs[seed_cID].find_element(
											cread_t(it_trg_db->rID));
						if (it_sd_db != _contigs[seed_cID].get_db_end()) {

							common_rID = it_trg_db->rID;
							to_inverse = (it_sd_db->is_fwd_strand ==
								it_trg_db->is_fwd_strand)? false : true;

							sd_dist_to_beg = it_sd_db->dist_to_beg;

							if (to_inverse) {
								trg_dist_to_beg = trg_clen -
										it_trg_db->dist_to_beg -
										it_trg_db->read_length;
							} else {
								trg_dist_to_beg = it_trg_db->dist_to_beg;
							}
						}
					}

				} else { /* new rID entry to be included */

					/* update rID_list for seed_cID */
					rID_list[seed_cID].insert(it_trg_db->rID);

					add_entries.push_back(*it_trg_db);
				}

				/* modify rID--> contigIdx mapping information */
				std::map<int, ipair_t>::iterator it_map =
						rID_cIDpair.find(it_trg_db->rID);
				if (it_map != rID_cIDpair.end()) {

					if (it_map->second.first == trg_cID){
						it_map->second.first = seed_cID;
						if (it_map->second.second == trg_cID) {
							it_map->second.second = seed_cID;
						}
					} else if (it_map->second.second == trg_cID){
						it_map->second.second = seed_cID;
					} else {
						abording("Contiger::merge_contig_cluster SC p1 failed");
					}
				}

			} // for (it_trg_db

			if (common_rID == -1) {
				abording ("Contiger::merge_contig_cluster SC p2 failed: "
						"no common rID b/t seed_cID and trg_cID");
			}
			/* now update add_entries according to [to_inverse],
			 * [sd_dist_to_beg] [trg_dist_to_beg] */
			int offset = sd_dist_to_beg - trg_dist_to_beg;
			for (int l = 0; l < (int) add_entries.size(); ++ l) {
				if (to_inverse) {
					add_entries[l].is_fwd_strand =
							add_entries[l].is_fwd_strand? false:true;
					add_entries[l].dist_to_beg =
						trg_clen - add_entries[l].dist_to_beg -
						add_entries[l].read_length + offset;

				} else add_entries[l].dist_to_beg += offset;

				min_dist = std::min(min_dist, offset);
			} // for (l

			all_entries.insert(all_entries.end(), add_entries.begin(),
					add_entries.end());

			/* destroy [trg_cID] */
			_contigs[trg_cID].clear();
			rID_list[trg_cID].clear();
		} // for (j

		if (all_entries.size()) {
			/* update all_entries */
			if (min_dist < 0) {

				/* add the [seed_cID] to all_entries */
				all_entries.insert(all_entries.end(),
						_contigs[seed_cID].get_db_begin(),
						_contigs[seed_cID].get_db_end());

				_contigs[seed_cID].clear();

				for (int j = 0; j < (int) all_entries.size(); ++ j) {
					all_entries[j].dist_to_beg -= min_dist;
				}
			}

			_contigs[seed_cID].insert(all_entries.begin(), all_entries.end());
			_contigs[seed_cID].calculate_length();
		}

	} // for (i
} // merge_contigs_cc

/*// OpenMP version of initializing contigs in batches.
void Contiger::batch_build_contig (const iivec_t& common_shingle_indices,
		int start, int end, const std::vector<shingle_t>& shingles,
		const std::vector<read_t>& reads) {

  #pragma omp parallel
  {
  std::vector<Contig> localContigs;
 	#pragma omp for

	for (int i = start; i <= end; ++ i) {

		// for each cluster, create vector of contig read format "cread_t"
		// as the initial cluster (done by candidate.init)


		Contig candidate;
		candidate.init (common_shingle_indices[i], shingles, reads);

		localContigs.push_back(candidate);

	} // for

	#pragma omp critical
	{
		_contigs.insert(_contigs.end(), localContigs.begin(), localContigs.end());
	}
  } //#pragma omp parallel
}
*/
