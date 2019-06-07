//============================================================================
// Project     : Diversifier
// Name        : Util.h
// Author      : Xiao Yang
// Created on  : Jul 14, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef UTIL_H_
#define UTIL_H_

// #include <sys/resource.h>  Getrusage()

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <omp.h>
#include <iomanip>
#include "jaz/fasta_file.hpp"

#if defined (_MSC_VER)
  #include "xny/getWinTime.hpp"
#else
  #include <sys/time.h>
#endif

typedef std::vector<int> ivec_t;
typedef std::vector<char> cvec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::vector<uvec_t> uuvec_t;
typedef std::vector<std::string> strvec_t;
typedef std::pair<int, int> ipair_t;
typedef std::pair<bool, bool> bpair_t;
typedef std::pair<uint32_t, uint32_t> upair_t;
typedef std::pair<uint8_t, uint8_t> u8pair_t;
typedef std::pair<std::string, std::string> strpair_t;
typedef std::map<int, int> imap_t;
typedef std::map<uint32_t, uint32_t> umap_t;
typedef std::set<int> iset_t;



/* Get all reads from fasta file and store in [output]
 * usage:
 *  strvec_t reads;
 *  readFasta(std::back_inserter(reads), faFileNm);
*/
template <typename oIter>
void readFasta (oIter output, const std::string& filenm) {
	bIO::FASTA_input fasta_sr(filenm);
	if (fasta_sr.operator bool() == false) {
		std::cout << "open " << filenm << " failed, it exists? \n";
		exit(1);
	}

	typedef bIO::FASTA_input::value_type value_type;

	while (++fasta_sr == true) {
		const value_type& v = *fasta_sr;
		std::string header = v.first;
		*(output++) = v.second;
	}
}

/*
 * used while making unionfind structure
 */
template <typename T>
T find(T elem, std::vector<T>& vec){

    if (elem == vec[elem])
        return elem;
    else {
        vec[elem] = find (vec[elem], vec);
        return vec[elem];
    }
}

/*used while making final clusters
*/
template <typename T>
T clsfind(T elem, const std::vector<T>& vec){
    while (elem != vec[elem]) elem = vec[elem];
    return elem;
}


/* @brief	If seq is low complexity return true
 * Currently if the following two criteria are satisfied:
 * 1) a single nt counts for over 50%
 * 2) there exists a polymer with length >= 9
 */
inline bool islowComplexity(const std::string& seq) {
	ivec_t cnts(5,0); // a c g t n
	char pre;
	// the size of a tandem repeat, current 1 nt
	int tmp_tandemSize = 1, max_tandemSize = 0;
	if (!seq.empty()) {
		pre = toupper(seq.at(0));
		switch (pre) {
		case 'A': ++cnts[0];
				  break;
		case 'C': ++cnts[1];
				  break;
		case 'G': ++cnts[2];
				  break;
		case 'T': ++cnts[3];
				  break;
		default:  ++cnts[4];
				  pre = 'N';
				  break;
		}
	}

	for (int i = 1; i < (int) seq.length(); ++ i) {
		switch (toupper(seq.at(i))){
		case 'A':
			++ cnts[0];
			if (pre == 'A') ++ tmp_tandemSize;
			else {
				if (tmp_tandemSize > max_tandemSize) {
					max_tandemSize = tmp_tandemSize;
				}
				tmp_tandemSize = 1;
				pre = 'A';
			}
			break;
		case 'C':
			++ cnts[1];
			if (pre == 'C') ++ tmp_tandemSize;
			else {
				if (tmp_tandemSize > max_tandemSize) {
					max_tandemSize = tmp_tandemSize;
				}
				tmp_tandemSize = 1;
				pre = 'C';
			}
			break;
		case 'G':
			++ cnts[2];
			if (pre == 'G') ++ tmp_tandemSize;
			else {
				if (tmp_tandemSize > max_tandemSize) {
					max_tandemSize = tmp_tandemSize;
				}
				tmp_tandemSize = 1;
				pre = 'G';
			}
			break;
		case 'T':
			++ cnts[3];
			if (pre == 'T') ++ tmp_tandemSize;
			else {
				if (tmp_tandemSize > max_tandemSize) {
					max_tandemSize = tmp_tandemSize;
				}
				tmp_tandemSize = 1;
				pre = 'T';
			}
			break;
		default:
			++ cnts[4];
			if (pre == 'N') ++ tmp_tandemSize;
			else {
				if (tmp_tandemSize > max_tandemSize) {
					max_tandemSize = tmp_tandemSize;
				}
				tmp_tandemSize = 1;
				pre = 'N';
			}
			break;
		}
	}
	max_tandemSize = std::max(tmp_tandemSize, max_tandemSize);

	int maxCnt = *std::max_element(cnts.begin(), cnts.end());
	if (100 * maxCnt/seq.length() >= 50 && max_tandemSize >= 9) {
		return true;
	}
	return false;
} // islowComplexity


/** Class that generates all possible B(n,k) combinations
   *  of the set [0,n). It requires O(k) memory.
   */
class Combination {
  public:
      /**
       */
      typedef unsigned int value_type;

      /** Create Combination object and set it to store first
       *  k-element combination over n-element set.
       */
      Combination(unsigned int n, unsigned int k) { reset(n, k); }

      /** Reset to first k-element combination over n-element set.
       */
      void reset(unsigned int n, unsigned int k) {
	  if (k > n/2) { std::cout << "you can specify different m to speed up n choose m\n"; }
	  n_ = n;
	  k_ = k;
	  k1_ = k - 1;
	  nk_ = n - k;
	  data_.resize(k_);
	  for (unsigned int i = 0; i < k_; ++i) data_[i] = i;
      } // reset

      /** Generate next combination (in the lexicographical order).
       *  @return true if next combination is available, false otherwise.
       */
      bool next() {
	  unsigned int i = k1_;
	  while ((i > 0) && (data_[i] == nk_ + i)) --i;
	  if ((i == 0) && (data_[i] == nk_)) return false;
	  ++data_[i];
	  for (; i < k1_; ++i) data_[i + 1] = data_[i] + 1;
	  return true;
      } // next

      /** @ return i-th element of the current combination.
       */
      value_type operator[](unsigned int i) const { return data_[i]; }

  private:
      unsigned int n_;
      unsigned int k_;
      unsigned int k1_;
      unsigned int nk_;
      std::vector<value_type> data_;

}; // class Combination

template <typename Iter>
int max_index (Iter begin, Iter end) {
	Iter max_iter = begin;
	for (Iter it = begin; it < end; ++ it) {
		if ((*it) > (*max_iter)) max_iter = it;
	}
	return (max_iter - begin);
}
/* create histogram of elements stored in [vec]
 * and write to [out] (&std::cout)
 */
inline void histogram (const std::vector<int>& vec, std::ostream* out){

	typedef std::map<int, int> imap_t;
	imap_t freq;
	for (int i = 0; i < (int) vec.size(); ++ i) {
		int entry = vec[i];
		imap_t::iterator it_map(freq.find(entry));
		if (it_map != freq.end()) {
			++it_map->second;
		} else {
			freq[entry] = 1;
		}

	}
	for (imap_t::iterator it_map = freq.begin();
			it_map != freq.end(); ++it_map) {
		(*out) << "\t" << it_map->first << "\t" << it_map->second << "\n";
	}
}

// some debugging print out scripts

inline void print_2dvecPairs (const std::vector<std::vector<ipair_t> >& twodvecP){
	for (int i = 0; i < (int) twodvecP.size(); ++ i) {
		std::cout << i << "\n";
		for (int j = 0; j < (int) twodvecP[i].size(); ++ j) {
			std::cout << "(" << twodvecP[i][j].first << ", "
					<<  twodvecP[i][j].second << ")\t";
		}
		std::cout << "\n";
	}
}

// print to file
inline void print_2dvec (const iivec_t& myVec, std::ostream* out){
	iivec_t::const_iterator it1 (myVec.begin());
	std::vector<int>::const_iterator it2;
	int index = 0;
	for (; it1 != myVec.end(); ++ it1){
		(*out) << (index ++) << ": ";
		for (it2 = it1->begin(); it2 != it1->end(); ++ it2){
			(*out) << (*it2) << "\t";
		}
		(*out) << "\n";
	}
}

inline void printHex(int i) {
	std::cout << " " << std::hex << i;
}

inline void print_1dvec (const ivec_t& myVec, std::ostream* out) {
	for (int i = 0; i < (int) myVec.size(); ++ i) {
		(*out) << "\t(" << i << "," << myVec[i] << ")";
		if ((i!=0) && (i % 6 == 0)) (*out) << "\n";
	}
	(*out) << "\n";
}

inline void abording (const std::string& msg) {
	std::cout << "[EXIT]: " << msg << "\n";
	exit(1);
}
inline void warning(const std::string& msg) {
	std::cout << "[WARNING]: " << msg << "\n";
}


/* Divide k (k <= 32) positions [0... k-1] into c chunks,
 * then mask (by setting corresponding bits to 0) all positions in
 * d chunks from each (c choose d) combination, store the results in masks
*/
template <typename int_t>
void masking (std::vector<int_t>& masks, int c, int k, int d) {

    if (c > k || c <= d || k > 32){
    		abording (" #chunks has to be (d, k] and k <= 32\n");
    }
    //1. randomize k kmer indices, and split into c chunks;

	// a) random shuffle indices
    ivec_t indices(k, 0);
    for (int i = 0; i < k; ++i) indices[i] = i;
    std::random_shuffle(indices.begin(), indices.end());

    //print out
    /*
    std::cout << "indices: \n";
    print_1dvec (indices, &std::cout);
	*/
    // b) segment indices and store each segment in chunks
    int unitSize = k/c; // unit size for each chunk
    if ( (float) k/c > k/c + 0.5) ++ unitSize;

    iivec_t chunks (c, ivec_t());
    for (int i = 0; i < c; ++ i) {
    		ivec_t::iterator start = indices.begin() + i * unitSize,
    				end = indices.begin() + (i + 1)* unitSize;
    		if (i == (c - 1)) end = indices.end();
    		chunks[i].insert(chunks[i].end(), start, end);
    }

    //print out
    /*
    std::cout << "chunks:\n";
    	print_2dvec (chunks, &std::cout);
	*/
    	// c) generate c choose d chunks, and for each combination, store
    	// all corresponding coordinates in maskIdx
    	int num = 0;
    	iivec_t maskIdx;

    	if (d > 0) {
		Combination cmb (c, d);
		do {
				ivec_t tmp;
				for (int i = 0; i < d; ++ i) {
					tmp.insert(tmp.end(), chunks[cmb[i]].begin(), chunks[cmb[i]].end());
				}
				maskIdx.push_back(tmp);
				++ num;
		} while (cmb.next() == true);

		// print out
		/*
		std::cout << "maskIdx:\n";
		print_2dvec (maskIdx, &std::cout);
		*/
		// 2. generate masks for each k-spectrum duplication
		masks.resize(num);
		for (int i = 0; i < num; masks[i++] = -1);

		for (int i = 0; i < num; ++ i) {
			for (int j = 0; j < (int) maskIdx[i].size(); ++j) {
				masks[i] &= ~(1 << (2 * maskIdx[i][j]));
				masks[i] &= ~(1 << (2 * maskIdx[i][j] + 1));
			}
		}
    	}
    	else {
    		masks.push_back(-1);
    	}
    //print out
    /*
    std::cout << "masks: \n";
    std::for_each(masks.begin(), masks.end(), printHex);
    std::cout << "\n";
	*/
}

inline double get_time() {
      timeval t;
      gettimeofday(&t, 0);
      return t.tv_sec + (0.000001 * t.tv_usec);
}

inline void print_time (const std::string& msg, double& timing){
    double cur_time = get_time();
    std::cout << "\n" << msg << " *** " << (cur_time - timing)/60 << " mins ***\n\n";
    timing = cur_time;
}

/*
 * Reading [batchSize] number of reads from files fH1 (and fH2 if specified)
 * [output] is the container type of [read_t]
 *
 * If fH2 is specified, then fH1 and fH2 are paired files and reading
 * are done alternatively.
 *
 * Assumptions: the starting position of fH1 and fH2 are the beginning
 * of each read record in the fasta or fastq format; and one sequence/line
 * no empty lines in the middle of file
 *
 *
 * Testing code
 * 	std::vector<read_t> reads;

	getBatchReads (std::back_inserter(reads), 10, myPara.inputFiles[0].handle,
			myPara.inputFiles[1].handle);

	for (int i = 0; i < (int) reads.size(); ++ i){
			std::cout << reads[i].header << "\n";
			std::cout << reads[i].read << "\n";
			std::cout << reads[i].score << "\n";
	}
	exit(1);

template <typename oIter>
void getBatchReads (oIter output, int batchSize, std::ifstream* fH1,
		std::ifstream* fH2) {

	std::string buf;
	read_t entry1, entry2;
	int count = 0;

	while (!fH1->eof()) {

		if (count > batchSize) return;

		if (!std::getline((*fH1), entry1.header).eof() &&
			!std::getline((*fH1), entry1.read).eof()) {

			count ++ ;

			if (fH2) { // paired file

				count ++ ;

				if (!std::getline((*fH2), entry2.header).eof() &&
					!std::getline((*fH2), entry2.read).eof()) {

					if (entry1.header.at(0) == '@'
							&& entry2.header.at(0) == '@') { //fastq
						// read two more lines
						if (std::getline((*fH1), buf).eof() ||
							std::getline((*fH1), entry1.score).eof() ||
							std::getline((*fH2), buf).eof() ||
							std::getline((*fH2), entry2.score).eof()) {
							std::cout << "getBatchReads: inconsistent read "
								<< "records in file pairs: (fh1, fh2) = "
								<< "(" << fH1 << "," << fH2 << ") ... exit\n";
							exit (1);
						}
						*(output++) = entry1;
						*(output++) = entry2;

					} else if (entry1.header.at(0) != '>' ||
							 entry2.header.at(0) != '>') {//fasta
						std::cout << "getBatchReads: inconsistent first "
							<<	"character of file pairs: (fh1, fh2) = "
							<< "(" << fH1 << "," << fH2 << ") ... exit\n";
						exit (1);
					}

				} else { // reading error
					std::cout << "getBatchReads: reading unaligned read "
						"record in second paired file (handle = " << fH2
						<< " )... exit\n";
					exit(1);
				}

			} else { // unpaired file

				if (entry1.header.at(0) == '@') { // fastq
					if (std::getline((*fH1), buf).eof() ||
						std::getline((*fH1), entry1.score).eof()) { // reading err
						std::cout << "getBatchReads: inconsistent read "
								<< "records in fastq file: fh = "
								<< fH1 << " ... exit\n";
						exit(1);
					}
					*(output++) = entry1;
				} else if (entry1.header.at(0) != '>') { // reading error
					std::cout << "getBatchReads: inconsistent first "
						<<	"character of file fh1 = "
						<< fH1 << " ... exit\n";
					exit (1);
				}
			}
		} else { // reading error

			// skip empty lines in the beginning and end of the file
			if (entry1.header.empty()) continue;

			std::cout << "getBatchReads: read record is unaligned "
						"in the first paired file (handle = " << fH1
						<< ")... exit\n";
			exit(1);
		}
	} // while
}

  extract all kmers (including duplicated ones) from [seq]
 from both forward and reverse complementary strands, only retain
 alphabetically smaller kmer as the representative, merge to [kmers] (unsorted)
 Test code:
	uvec_t kArray;
	std::string str = "TTTTACGTtccat";
	extractKmerSequence (kArray, str, 4);
	std::cout << str << "\n\n";
	for (int i = 0; i < kArray.size(); ++ i) {
		std::cout << str.substr(i, 4) << "\n";
		std::cout << toString (kArray[i], 4) << "\n\n";
	}
	exit(1);

template <typename container_t>
void KmerExtr (container_t& kmers, const std::string& seq, int k) {

	typedef typename container_t::value_type elem_t;

	// lowest 2k bits = 1
    elem_t MaskLowerKbits = pow(2, 2*k) - 1;
    container_t MaskTop2bits(4, 0);
    for (elem_t i = 1; i < 4; ++ i){
        MaskTop2bits[i] = (i << (2*(k - 1)));
    }
    std::reverse(MaskTop2bits.begin(), MaskTop2bits.end()); // t g c a

    char* addr = const_cast<char*> (seq.c_str());
    int len = (int) seq.length() - k + 1;

	elem_t ID = 0, rcID = 0;
	int j = 0;
	bool flag = false; //true: the previous kmer contains no Ns
	while (j < len){
		if (flag){
			int c = char2bits (*(addr + j + k -1));
			if (c == -1) {
				j += k;
				flag = false;
			}
			else {
				ID  = ((ID << 2 | c) & MaskLowerKbits);
				rcID = (rcID >> 2);
				rcID |= MaskTop2bits[c];

				if (rcID < ID)
					kmers.push_back(rcID);
				else kmers.push_back(ID);
				++ j;
			}
		}
		else{
			if (toID (ID, addr + j, k)) {
				rcID = reverseComplBits <elem_t, int> (ID, k);
				flag = true;
				if (rcID < ID)
					kmers.push_back(rcID);
				else kmers.push_back(ID);
			}
			++ j;
		}
	}
}

 extract (kmer, pos) tuple from seq, including both fwd and rvc strands
 of seq; for rvc strand, the kmer pos is w.r.t to the fwd strand but with
 a negative sign, e.g. seq = acggttc, then "acc" occurs at position -2

Test code:
	std::vector<std::pair<uint32_t, int> > kArray;
	std::string str = "TTTTACGTtccat";
	KPExtr (kArray, str, 4);
	std::cout << str << "\n\n";
	std::cout << "kArray size = " << kArray.size() << "\n";
	for (int i = 0; i < kArray.size()/2; ++ i) {
		std::cout << i << ":\n";
		std::cout << str.substr(i, 4) << "\n";
		std::cout << toString <uint32_t>(kArray[2*i].first, 4) << "\n";
		std::cout << toString <uint32_t>(kArray[2*i + 1].first, 4) << "\n\n";
	}
	exit(1);

template <typename int_t>
void KPExtr (std::vector<std::pair<int_t, int_t> >& kmers,
		const std::string& seq, int k) {

	// lowest 2k bits = 1
    int_t MaskLowerKbits = pow(2, 2*k) - 1;
    std::vector<int_t> MaskTop2bits(4, 0);
    for (int_t i = 1; i < 4; ++ i){
        MaskTop2bits[i] = (i << (2*(k - 1)));
    }
    std::reverse(MaskTop2bits.begin(), MaskTop2bits.end()); // t g c a

    char* addr = const_cast<char*> (seq.c_str());
    int len = (int) seq.length() - k + 1;

    int_t ID = 0, rcID = 0;
	int j = 0;
	bool flag = false; //true: the previous kmer contains no Ns
	while (j < len){
		if (flag){
			int c = char2bits (*(addr + j + k -1));
			if (c == -1) {
				j += k;
				flag = false;
			}
			else {
				ID  = ((ID << 2 | c) & MaskLowerKbits);
				rcID = (rcID >> 2);
				rcID |= MaskTop2bits[c];

				kmers.push_back(std::pair<int_t, int_t> (ID, j));
				kmers.push_back(std::pair<int_t, int_t> (rcID, (-1)*j));

				++ j;
			}
		}
		else{
			if (toID (ID, addr + j, k)) {
				rcID = reverseComplBits <int_t, int_t> (ID, k);
				flag = true;

				kmers.push_back(std::pair<int_t, int_t> (ID, j));
				kmers.push_back(std::pair<int_t, int_t> (rcID, (-1)*j));
			}
			++ j;
		}
	}
}*/
#endif /* UTIL_H_ */
