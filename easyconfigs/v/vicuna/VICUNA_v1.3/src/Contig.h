//============================================================================
// Project     : Diversifier
// Name        : Contig.h
// Author      : Xiao Yang
// Created on  : Aug 26, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef CONTIG_H_
#define CONTIG_H_

#include "xutil.h"
#include "ReadIndexer.h"
#include "xny/seq_manip.hpp"
#include "xny/seq_pair_manip.hpp"

struct shingle_t{
	int rID;
	uint32_t fwd_ss;  // fwd_super_shingle
	uint32_t rv_ss;
	uint8_t fwd_min_pos; // fwd_min_shingle position w.r.t fwd strand
	uint8_t rv_min_pos;  // rv_min_shingle position w.r.t rvc strand
	bool is_paired;
};

/* Contig Read Type */
struct cread_t {
	cread_t (){ rID = -1; }
	cread_t (int id): rID(id) {}
	uint8_t indel_length () const {
		uint8_t indel_len = 0;
		for (int i = 0; i < (int) indels.size(); ++ i){
			indel_len += indels[i].second;
		}
		return indel_len;
	}
	int rID;			    // read ID
	int dist_to_beg; 	// read start position w.r.t. contig start
	uint8_t read_length;
	bool is_fwd_strand; // w.r.t. the input strand
	bool is_paired;
	std::vector<u8pair_t> indels; // indel: <start, len>
};

/* compare cread_t by rID [default] */
inline bool operator<(const cread_t& lhs, const cread_t& rhs) {
	return lhs.rID < rhs.rID;
}

/* sort ipair_t in increasing order by the first element */
struct comp_ipair {
	bool operator () (const ipair_t& lhs, const ipair_t& rhs) const {
		return lhs.first < rhs.first;
	}
};

/* sort u8pair_t in increasing order by the first element */
struct cmp_u8pair {
	bool operator () (const u8pair_t& lhs, const u8pair_t& rhs) const {
		return lhs.first < rhs.first;
	}
};

class Contig {
public:
	Contig () { length_ = 0; };

	virtual ~Contig() {};
	inline int number_of_reads () const { return db_.size(); }
	inline void set_db (const std::set<cread_t>& db) { db_ = db; }
	inline std::set<cread_t>::iterator get_db_begin ()
											{ return db_.begin(); }
	inline std::set<cread_t>::iterator get_db_end ()
											{ return db_.end(); }
	inline std::set<cread_t>::const_iterator get_db_begin ()
										const { return db_.begin(); }
	inline std::set<cread_t>::const_reverse_iterator get_db_rbegin ()
										const { return db_.rbegin(); }
	inline std::set<cread_t>::const_iterator get_db_end ()
										const { return db_.end(); }
	inline std::set<cread_t>::const_iterator
		find_element(const cread_t& read) const { return db_.find(read); }
	inline std::set<cread_t>::iterator
		find_element(const cread_t& read) { return db_.find(read); }

	inline void set_length(int length) { length_ = length; }
	inline int get_length () const { return length_; }
	void reset_db ();
	void update_dist_to_beg (int dist) {
		if (dist == 0) return;
		std::set<cread_t> tmp_db;
		std::set<cread_t>::iterator it = db_.begin();
		for (; it != db_.end(); ++ it) {
			cread_t tmp_read = (*it);
			tmp_read.dist_to_beg += dist;
			tmp_db.insert(tmp_read);
		}
		db_ = tmp_db;
	} // update_dist_to_beg

	void sort_by_dist (std::vector<ipair_t>& dist_to_rID) const;

	template <typename iIter>
	inline void insert (iIter begin, iIter end) {
		for (iIter it = begin; it != end; ++ it) db_.insert(*it);
	}
	inline void insert (const cread_t& cread) { db_.insert(cread); }
	inline void clear () { db_.clear(); length_ = 0; }

	void calculate_length () {
		length_ = 0;
		std::set<cread_t>::iterator it = db_.begin();
		for (; it != db_.end(); ++ it) {
			length_ = std::max (length_, int (it->read_length +
					it->dist_to_beg + it->indel_length()));
		}
	}

	void get_rIDs (iset_t& rIDs) {
		std::set<cread_t>::iterator it = db_.begin();
		for (; it != db_.end(); ++ it) rIDs.insert(it->rID);
	}
	void get_rIDs (ivec_t& rIDs) {
		std::set<cread_t>::iterator it = db_.begin();
		for (; it != db_.end(); ++ it) rIDs.push_back(it->rID);
	}

	void init (const ivec_t& shingle_indices,
			   const std::vector<shingle_t>& shingles,
			   const std::vector<read_t>& reads);

	Contig split (const std::set<ipair_t, comp_ipair>& rhs_indices);

	bool validate (const std::set<ipair_t, comp_ipair>& rIndices,
		const std::vector<read_t>& reads, std::set<ipair_t, comp_ipair>& diverged,
		int overhang, int divergence);


	bool merge_via_rID (Contig& rhs, int common_rID);
	void merge_via_overlap (const xny::align_t& algn, Contig& rhs);

	void print_alignment (std::ofstream& oHandle, int contigID,
			const ReadIndexer& rIndexer, const Trimmer& myTrimmer);
	void print (int contigID, const ReadIndexer& rIndexer, const Trimmer& myTrimmer,
			bool to_print_profile, bool to_print_reads, int tail_length = -1);
	void print_profile(const std::string& consensus, const ivec_t& profile) const;
	void print_reads (const std::vector<read_t>& reads, int tail_length) const;

	/* normal profile */
	void generate_profile (std::string& consensus, ivec_t& profile,
			const std::set<ipair_t, comp_ipair>& rIndices,
			const std::vector<read_t>& reads) const;
	void update_profile (std::string& consensus, ivec_t& profile,
		 std::set<cread_t>::const_iterator it, const std::string& read_seq,
		 int value) const;

	/* weighted profile based on distance to read beginning */
	void generate_weighted_profile (std::string& consensus,
		ivec_t& profile, const std::set<ipair_t, comp_ipair>& rIndices,
		const std::vector<read_t>& reads) const;
	void update_weighted_profile (std::string& consensus, ivec_t& profile,
		 std::set<cread_t>::const_iterator it, const std::string& read_seq,
		 int sign) const;

	void inverse_all_reads ();

private:
	std::set<cread_t> db_;
	int length_;

	inline void update_dist_to_beg (std::vector<cread_t>& lhs, int dist) {
		if (dist == 0) return;
		for (int i = 0; i < (int) lhs.size(); ++ i) lhs[i].dist_to_beg += dist;
	}

	void update_all_reads (int extr_dist,
			const std::vector<ipair_t>& contig_gaps);
	void update_read (std::set<cread_t>::iterator it, int dist,
		const std::vector<u8pair_t>& gaps);

	bool verify_consensus (std::string& consensus, ivec_t& profile,
			const std::set<ipair_t, comp_ipair>& rIndices,
			const std::vector<read_t>& reads,
			std::set<ipair_t, comp_ipair>& diverged_rIndices,
			int overhang, int divergence);

	void update_profile_column (std::string& consensus, ivec_t& profile,
			int col, char newchar, int value) const;

	void inverse_read (int contig_length, cread_t& cread);
}; // class contig

struct comp_contig_ge {
	bool operator () (const Contig& lhs, const Contig& rhs) const {
		return rhs.get_length() < lhs.get_length();
	}
};

struct comp_contig_lt {
	bool operator () (const Contig& lhs, const Contig& rhs) const {
		return lhs.get_length() < rhs.get_length();
	}
};
#endif /* CONTIG_H_ */
