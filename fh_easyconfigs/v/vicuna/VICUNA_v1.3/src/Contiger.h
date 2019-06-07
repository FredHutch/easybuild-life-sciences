//============================================================================
// Project     : Diversifier
// Name        : Contiger.h
// Author      : Xiao Yang
// Created on  : Aug 29, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef CONTIGER_H_
#define CONTIGER_H_

#include "Parameter.h"
#include "xutil.h"
#include "Trimmer.h"
#include "ReadIndexer.h"
#include "Contig.h"
#include "Delegate.h"
#include "xny/seq_pair_manip.hpp"
#include "jaz/functional_add.hpp"
#include "jaz/string_add.hpp"
/*
 *
 */
/* neighbor of (this) contig */
struct contig_nb_t {
	contig_nb_t (int id): cID(id) { pp = pm = mp = mm = 0; }
	inline int sum () const { return pp + pm + mp + mm; }
	int cID;
	int pp;
	int pm;
	int mp;
	int mm;
} ;

struct cmp_contig_nb_lt {
	bool operator () (const contig_nb_t& lhs, const contig_nb_t& rhs)
		const {	return lhs.sum() < rhs.sum();
	}
};
struct cmp_contig_nb_id {
	bool operator () (const contig_nb_t& lhs, const contig_nb_t& rhs)
		const {	return lhs.cID < rhs.cID;
	}
};

/* comparing sh in function generate_ss */
struct comp_sh {
	bool operator () (const std::pair<uint32_t, uint8_t>& lhs,
		const std::pair<uint32_t, uint8_t>& rhs) const {
		return lhs.first < rhs.first;
	}
};

/* spaced-seed --> readID, pos, dir */
typedef struct SpacedSeed {
	SpacedSeed (){}
	SpacedSeed (int id, uint8_t p, bool d):
		rID (id), pos (p), dir (d) {}
	int rID;
	uint8_t pos;
	bool dir;
} spaced_seed_t;

class Contiger {
public:
	Contiger (const Parameter& myPara){
		_w1 = myPara.c_w1;
		_w2 = myPara.c_w2;
		_batchsz = myPara.batchSize;
		read_overhang_ = myPara.c_read_overhang;
		percent_diverge_ = myPara.c_percent_diverge;
		min_profile_col_weight_ = myPara.c_min_profile_col_weight;
		min_consensus_base_ratio_ = myPara.c_min_consensus_base_ratio;
		max_contig_overhang_ = myPara.c_max_contig_overhang;
		seed_kmer_len_ = myPara.c_seed_kmer_len;
    		min_contig_overlap_ = myPara.c_min_contig_overlap;
    		_min_contig_links = myPara.c_min_contig_links;
    		min_identity_ = myPara.c_min_identity;
    		min_perc_pol_ = myPara.c_min_perc_pol;
    		max_variant_len_ = myPara.c_max_variant_len;
    		_min_output_contig_len = myPara.min_output_contig_len;
    		_out_dir = myPara.oDIRNm;
	}
	virtual ~Contiger(){};
	void run (bool is_filter_on, const ivec_t& binned_rIDs,
			const Trimmer& myTrimmer, const ReadIndexer& rIndexer);
	void scaffold ();
	inline int num_contigs () const { return _contigs.size(); }
	//inline void add_contig (const Contig& rhs) { _contigs.push_back(rhs);}
	inline int number_of_reads (int ID) {
								return _contigs[ID].number_of_reads(); }
private:

	std::vector<Contig> _contigs;
	/* input parameters */
	int _w1, _w2, _batchsz, read_overhang_, percent_diverge_,
		min_profile_col_weight_, min_consensus_base_ratio_,
		max_contig_overhang_, seed_kmer_len_, min_contig_overlap_,
		_min_contig_links, min_identity_, min_perc_pol_, max_variant_len_,
		_min_output_contig_len;
	std::string _out_dir;

	/* construct super shingle structure */
	void super_shingling (std::vector<shingle_t>& shingles,
		const ivec_t& rIDs, const ReadIndexer& rIndexer,
		const Trimmer& myTrimmer);
	void batch_ss (std::vector<shingle_t>& shingles,
			const std::vector<read_t>& reads, const Trimmer& myTrimmer);
	bool generate_ss(uint32_t& ss, uint8_t& common_shingle_pos,
			const std::string& read);

	/* rely on super-shingles to cluster reads */
	void ss_clustering (iivec_t& ss2rIDs,
						const std::vector<shingle_t>& shingles);

	/* construct contigs */
	void init_contig (const iivec_t& common_shingle_indices,
		const std::vector<shingle_t>& shingles, const Trimmer& myTrimmer,
		const ReadIndexer& rIndexer);
	//void batch_build_contig (const iivec_t& common_shingle_indices,
	//		int beg, int end, const std::vector<shingle_t>& shingles,
	//		const std::vector<read_t>& reads);

	/* fishing reads and add to initial contigs via kmers */
	void fishing (bool is_filter_on, const Trimmer& myTrimmer,
			const ReadIndexer& rIndexer);
	void identify_reads_involved_in_fishing (ivec_t& rIDs,
			int input_read_cnt,	bool filtered);
	void batch_fishing (std::map<uint64_t, spaced_seed_t>& STable,
		imap_t& rID_cID_map, const std::vector<bool>& seed_template,
		const std::vector<read_t>& reads);
	void cluster_read_via_seed (cread_t& query_cread,
	 imap_t& rID_cID_map, int seed_len, const spaced_seed_t& query_entry,
		std::map<uint64_t, spaced_seed_t>::iterator it_hit_STable);

	void generate_seed (std::vector<bool>& seed);

	/* merge contigs using rID->contig_index structure, union find*/
	void merge_contig_via_rID_uf ();
	/* Alternatively, based on connected component */
	void merge_contig_via_rID_multi();
	void merge_contigs_cc_multi (std::map<int, ivec_t>& rID_cIDlist_map,
		std::vector<iset_t>& rID_list, const iivec_t& cluster_of_cIDs);
	void merge_contig_via_rID ();
	void merge_contigs_cc (std::map<int, ipair_t>& rID_cIDpair,
		std::vector<iset_t>& rID_list, const iivec_t& cluster_of_cIDs);
	/* merge contigs via overlaps */
	/*
	bool merged (int cID0, int cID1, bool is_rvc1, delegate_t& dg0,
		delegate_t& dg1, const std::string& rl0, const std::string& rl1);
	*/
	//void get_reliable_regions (strpair_t& rls, std::pair<std::vector<delegate_t>::iterator,
	//		std::vector<delegate_t>::iterator> it_dg, bool rvc1, int refcID,
	//		const contig_nb_t& nb, const std::vector<delegate_t>& delegates);
	bool merge_contig_via_overlap (int cID0, int cID1, bool rvc1,
		std::vector<Delegate>& delegates);
	void low_frq_var_rec (std::string& alignment, int& aln_start,
		int& aln_end, std::vector<ipair_t>::const_iterator it_lfv_start,
		std::vector<ipair_t>::const_iterator it_lfv_end, char ins);
	void close_alt_indels (xny::align_t& aln, const std::string& s0,
			const std::string& s1);

	void compute_add_nbs (std::vector<contig_nb_t>& add_nbs,
				const std::vector<contig_nb_t>& nbs,
				const std::vector<Delegate>& dgs);

	bool rvc_checker (const contig_nb_t& nb);
	bool is_in_nbs (int& idx, std::vector<contig_nb_t>& nbs, int cID);
	void merge_nbs (std::vector<contig_nb_t>& lhs_nbs,
					std::vector<contig_nb_t>& rhs_nbs);
	void update_rcMap (imap_t& rcMap, int src_cID, int tgt_cID);

	/* validating contigs */
	void validate_contig (const Trimmer& myTrimmer,
						  const ReadIndexer& rIndexer);
	void batch_validate_contig (int start, int end,
			const imap_t& rID_index, const std::vector<read_t>& reads,
			const Trimmer& myTrimmer,  const ReadIndexer& rIndexer);
	void divide_contig_by_n (int index, std::vector<Contig>& contigs);
	void divide_contig_by_minoverlap (int index,
				std::vector<Contig>& contigs, const imap_t& rID_index,
				const std::vector<read_t>& reads,
				const Trimmer& myTrimmer,
											   const ReadIndexer& rIndexer);
	void divide_contig_by_minoverlap_0 (int index,
				std::vector<Contig>& contigs, const imap_t& rID_index,
				const std::vector<read_t>& reads);
	/* sort & remove empty entries */
	void reorganize_contigs ();
	/* extend existing contigs using paired read information */
	void extend_contig (const Trimmer& myTrimmer,
					    const ReadIndexer& rIndexer);
	bool next_cID(int& start_idx);
	void generate_rID_cID_map (imap_t& rID_cID_map, bool paired_read_only);

	void get_contig_neighbors (std::vector<contig_nb_t>& nbs, int cID,
			const imap_t& rID_cID);
	void compute_delegates (std::vector<Delegate>& delegates,
			const std::vector<contig_nb_t>& nbs, const Trimmer& myTrimmer,
			const ReadIndexer& rIndexer);
	//void compute_delegate_i (delegate_t& delegate, const Contig& contig,
	//		const std::set<ipair_t, comp_ipair>& rIndices,
	//		const std::vector<read_t>& reads);
	void output (const std::vector<Delegate>& dgs,
			const Trimmer& myTrimmer, const ReadIndexer& rIndexer);

};


#endif /* CONTIGER_H_ */
