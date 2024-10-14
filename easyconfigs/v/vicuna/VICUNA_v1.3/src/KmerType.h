//============================================================================
// Project     : GenomeAnalysis
// Name        : KmerType.h
// Author      : Xiao Yang
// Created on  : Jul 25, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef KMERTYPE_H_
#define KMERTYPE_H_

#include <inttypes.h>

/*
 * usage: g++ -DLONG_KMERS test.cpp -o test.exe
 */
#ifdef LONGKMER
#define kmer_t uint64_t
#else
#define kmer_t uint32_t
#endif


#endif /* KMERTYPE_H_ */
