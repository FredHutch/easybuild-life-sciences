//============================================================================
// Project     : Diversifier
// Name        : Delegate.h
// Author      : Xiao Yang
// Created on  : Oct 26, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef DELEGATE_H_
#define DELEGATE_H_

#include "xutil.h"
#include "xny/seq_manip.hpp"
#include "xny/seq_pair_manip.hpp"

/* The class to retain "consensus", "profile", and reliable interval
 */

class Delegate {
public:
	Delegate () { cID = -1; }
	Delegate (int id): cID(id) {}
	inline void set_consensus (const std::string& cons) { consensus = cons;}
	inline void set_profile (const ivec_t& prfl) { profile = prfl;}
	inline void set_cID (int id) { cID = id;}
	inline int get_cID () const { return cID; }
	inline int lfv_size () const { return lfv.size(); }
	inline int length () const { return consensus.length(); }
	inline std::string get_consensus() const { return consensus; }
	inline std::vector<ipair_t>::const_iterator get_lfv_begin ()
												{ return lfv.begin(); }
	inline std::vector<ipair_t>::const_iterator get_lfv_end ()
												{ return lfv.end(); }
	void rvc ();
	void merge (const xny::align_t& aln, const Delegate& rhs);
	void id_low_frq_var (int min_perc_pol, int max_variant_len);
	void low_frq_var_free_consn (std::string& rc) const;
	void print_profile();
private:
	int cID;
	std::string consensus;
	ivec_t profile;
	std::vector<ipair_t> lfv; /* low freq variants <pos, len>,
	 	 	 	 	 	 	 	 where pos is wrt consensus */

	/* the following component is implemented and can be turned on */
	//ipair_t ri; /* reliable interval of the consensus */
	//void get_reliable_region(int min_prfl_col_weight,
	//	int min_cons_base_ratio, int max_contig_overhang);
};

struct cmp_delegate {
	bool operator () (const Delegate& lhs, const Delegate& rhs) const {
		return lhs.get_cID() < rhs.get_cID();
	}
};

#endif /* DELEGATE_H_ */
