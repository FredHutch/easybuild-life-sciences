//============================================================================
// Project     : DivAnalysis
// Name        : readEntry.h
// Author      : Xiao Yang
// Created on  : Dec 15, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef READENTRY_H_
#define READENTRY_H_

#include <iostream>
#include <string>

struct rentry_t {
	rentry_t(int id, const std::string& h, std::string r, std::string s):
		ID (id), header(h), read(r), score(s){}
	rentry_t(int id): ID(id) {}
	rentry_t(){}
	bool paired;
	int ID;
	std::string header;
	std::string read;
	std::string score;
};


#endif /* READENTRY_H_ */
