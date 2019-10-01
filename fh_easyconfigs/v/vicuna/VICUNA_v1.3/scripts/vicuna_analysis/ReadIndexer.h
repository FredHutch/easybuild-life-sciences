//============================================================================
// Project     : DivAnalysis
// Name        : ReadIndexer.h
// Author      : Xiao Yang
// Created on  : Aug 18, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================

#ifndef READINDEXER_H_
#define READINDEXER_H_

#include <vector>
#include <string>
#include <inttypes.h>

#include "Parameter.h"
#include "readEntry.h"
#include "Trimmer.h"
#include "xny/file_manip.hpp"


typedef struct ifiles {
	ifiles(const std::string& fn, uint32_t num, bool f1, bool f2) :
			filename(fn), numReads(num), fastq(f1), paired(f2) {
	}
	std::string filename;
	std::ifstream* handle;
	uint32_t numReads;
	bool fastq;
	bool paired;
} file_t;

struct ReadIdx {
	ReadIdx(int fid, uint64_t os) :
			fileID(fid), offset(os) {
	}
	int fileID;
	uint64_t offset;
};


class ReadIndexer {
public:
	ReadIndexer(const std::string& pFqDir, const std::string& npFqDir,
			const std::string& pFaDir, const std::string& npFaDir) {

		getInputFilelist(pFqDir, npFqDir, pFaDir, npFaDir);
		create();
		int totalReads = 0;
		for (int i = 0; i < (int) filelist_.size(); ++ i) {
			totalReads += filelist_[i].numReads;
			std::cout << "\t"
					  << xny::get_suffix(filelist_[i].filename, true, "/")
					  << "\t#reads: " << filelist_[i].numReads << "\n";
		}
		std::cout << "\tTotal # Input Reads: " << totalReads << "\n\n";
	}
	~ReadIndexer() {
		for (int i = 0; i < (int) filelist_.size(); ++i) {
			filelist_[i].handle->close();
			delete (filelist_[i].handle);
		}
	}

	template<typename oIter, typename iIter>
	void get_trimmed_reads (oIter output, iIter first, iIter last,
			const Trimmer& myTrimmer) const {
		// get all involved fileIDs
		for (; first != last; ++first) {

			if ((*first) >= (int) _rIdx.size()) {
				std::cout << "\t[warning] getReads: " << (*first)
						<< " is out of range ... stop retrieving\n";
				return;
			}
			rentry_t entry;
			int fID = _rIdx[*first].fileID;
			entry.ID = (*first);
			entry.paired = isPaired(*first);
			seekRead (entry, _rIdx[(*first)].offset, filelist_[fID].handle);
			myTrimmer.trim (entry.read, entry);
			(*output) = entry;
			++ output;
		}
	} // get trimmed reads

	bool isPaired (int rID) const {
		if (rID >= (int) _rIdx.size()) {
			std::cout << "[Warning] isPaired: rID " << rID
					<< " is out of range\n";
			return false;
		}
		int fID = _rIdx[rID].fileID;
		if (fID >= (int) size()) {
			std::cout << "isPaired: Sanity check failed ...exit\n";
			exit(1);
		}
		return filelist_[fID].paired;
	} // isPaired

	uint64_t size() const {
		return (_rIdx.size());
	}
private:
	std::vector<ReadIdx> _rIdx;
	std::string _buf, _buf_rv;
	uint64_t _is_pos, _is_pos_rv;
	std::string pFqDir, npFqDir, pFaDir, npFaDir;
	std::vector<file_t> filelist_;

	void seekRead(rentry_t& entry, uint64_t offset, std::ifstream* fH)
	const {
		std::string dummie;
		if (!fH->good()) fH->clear();
		fH->seekg(offset, std::ios_base::beg);
		std::getline((*fH), entry.header);
		std::getline((*fH), entry.read);
		if (!entry.header.empty()) {
			if (entry.header.at(0) == '@') {
				std::getline((*fH), dummie); // dummie variable
				std::getline((*fH), entry.score);
			} else if (entry.header.at(0) != '>') {
				std::cout << "seekRead problem: line starting with"
						<< entry.header.at(0) << " at offset " << offset
						<< " exit ! \n";
				exit(1);
			}
		} else {
			std::cout << "seekRead problem: line is empty, exit! \n";
			exit(1);
		}
	} //seekRead

	void getInputFilelist(const std::string& pFqDir,
			const std::string& npFqDir, const std::string& pFaDir,
			const std::string& npFaDir) {

		if (pFqDir.empty() && npFqDir.empty() && pFaDir.empty()
				&& npFaDir.empty()) {
			std::cout << "\n\nError: no input directory specified !!!\n";
			exit(1);
		} else {
			std::vector<std::string> fileNms;
			if (!pFqDir.empty()) {
				if (xny::dir_list(pFqDir, std::back_inserter(fileNms))
															== false) {
					std::cout << "can't read from dir: " << pFqDir << "\n";
					exit(1);
				} else {
					for (int i = 0; i < (int) fileNms.size(); ++i) {

						if (fileNms[i].compare(".") == 0
								|| fileNms[i].compare("..") == 0)
							continue;

						if (xny::get_suffix(
								fileNms[i], true, ".").compare("fq") == 0
						 ||
						 xny::get_suffix(fileNms[i],
								 true, ".").compare("fastq") == 0){

							filelist_.push_back(
								file_t(pFqDir + "/" + fileNms[i], 0, 1, 1));
							filelist_.back().handle = new std::ifstream(
									filelist_.back().filename.c_str());
							if (!filelist_.back().handle->good()) {
								std::cout << "Can't open file "
										<< filelist_.back().filename << "\n";
								exit(1);
							}
						}
					}
				}
			} // if
			fileNms.clear();
			if (!pFaDir.empty()) {
				if (xny::dir_list(pFaDir, std::back_inserter(fileNms))
															== false) {
					std::cout << "can't read from dir: " << pFaDir << "\n";
					exit(1);
				} else {
					for (int i = 0; i < (int) fileNms.size(); ++i) {
						if (fileNms[i].compare(".") == 0
								|| fileNms[i].compare("..") == 0)
							continue;
						if (xny::get_suffix(fileNms[i],
								true, ".").compare("fa") == 0
									||
							xny::get_suffix(fileNms[i],
								 true, ".").compare("fasta") == 0){

							filelist_.push_back(
								file_t(pFaDir + "/" + fileNms[i], 0, 0, 1));
							filelist_.back().handle = new std::ifstream(
									filelist_.back().filename.c_str());
							if (!filelist_.back().handle->good()) {
								std::cout << "Can't open file "
										<< filelist_.back().filename << "\n";
								exit(1);
							}
						}
					}
				}
			} // if
			fileNms.clear();
			if (!npFqDir.empty()) {
				if (xny::dir_list(npFqDir, std::back_inserter(fileNms)) == false) {
					std::cout << "can't read from dir: " << npFqDir << "\n";
					exit(1);
				} else {
					for (int i = 0; i < (int) fileNms.size(); ++i) {

						if (fileNms[i].compare(".") == 0
								|| fileNms[i].compare("..") == 0)
							continue;

						if (xny::get_suffix(fileNms[i],
								true, ".").compare("fq") == 0
								||
							xny::get_suffix(fileNms[i],
								 true, ".").compare("fastq") == 0){

							filelist_.push_back(
								file_t(npFqDir + "/" + fileNms[i],
										0, 1, 0));
							filelist_.back().handle = new std::ifstream(
									filelist_.back().filename.c_str());
							if (!filelist_.back().handle->good()) {
								std::cout << "Can't open file "
									<< filelist_.back().filename << "\n";
								exit(1);
							}
						}
					}
				}
			} // if
			fileNms.clear();
			if (!npFaDir.empty()) {
				if (xny::dir_list(npFaDir, std::back_inserter(fileNms))
															== false) {
					std::cout << "can't read from dir: "
							  << npFaDir << "\n";
					exit(1);
				} else {
					for (int i = 0; i < (int) fileNms.size(); ++i) {

						if (fileNms[i].compare(".") == 0
								|| fileNms[i].compare("..") == 0)
							continue;

						if (xny::get_suffix(fileNms[i],
								true, ".").compare("fa") == 0
								||
							xny::get_suffix(fileNms[i],
									true, ".").compare("fasta") == 0){

							filelist_.push_back(
								file_t(npFaDir + "/" + fileNms[i],
										0, 0, 0));
							filelist_.back().handle = new std::ifstream(
									filelist_.back().filename.c_str());
							if (!filelist_.back().handle->good()) {
								std::cout << "Can't open file "
									<< filelist_.back().filename << "\n";
								exit(1);
							}
						}
					}
				}
			} // if
		} // else
	} //getInputFilelist

	void create() {

		std::vector<ReadIdx>().swap(_rIdx);

		for (int i = 0; i < (int) filelist_.size(); ++i) {

			std::cout << "Indexing ...\n\n";
			std::cout << "\t" << filelist_[i].filename << "\n";

			std::ifstream* fH = filelist_[i].handle;

			int numReads = 0;

			if (filelist_[i].paired) { // paired read files
				if (i + 1 < (int) filelist_.size()) {
					std::cout << "\t" << filelist_[i + 1].filename << "\n";
				} else {
					std::cout << "create err: No " << i + 1
							<< " th file, exit\n";
					exit(1);
				}
				std::ifstream* fH_rv = filelist_[i + 1].handle;

				while (!std::getline((*fH), _buf).eof()
						&& !std::getline((*fH_rv), _buf_rv).eof()) {

					numReads ++;

					if (!_buf.empty() && !_buf_rv.empty()) {
						if (_buf.at(0) == '@' && _buf_rv.at(0) == '@') { //fastq
							_is_pos = fH->tellg();
							_is_pos_rv = fH_rv->tellg();
							// getline trunks the last char
							_is_pos -= (_buf.size() + 1);
							_is_pos_rv -= (_buf_rv.size() + 1);
							_rIdx.push_back(ReadIdx(i, _is_pos));
							_rIdx.push_back(ReadIdx(i + 1, _is_pos_rv));
							for (int j = 0; j < 3; ++j) {
								std::getline((*fH), _buf);
								std::getline((*fH_rv), _buf_rv);
								if (_buf.empty() || _buf_rv.empty()){
									std::cout<< "create: err reading paired fastq files:\n"
											<< filelist_[i].filename << "\n"
											<< filelist_[i + 1].filename
											<< " ... exit \n";
									exit(1);
								}
							}
						} else if (_buf.at(0) == '>' && _buf_rv.at(0) == '>') {	//fasta
							_is_pos = fH->tellg();
							_is_pos_rv = fH_rv->tellg();
							_is_pos -= (_buf.size() + 1);
							_is_pos_rv -= (_buf_rv.size() + 1);
							_rIdx.push_back(ReadIdx(i, _is_pos));
							_rIdx.push_back(ReadIdx(i + 1, _is_pos_rv));

							std::cout << _buf << "\n";
							std::cout << _buf_rv << "\n\n";

							std::getline((*fH), _buf);
							std::getline((*fH_rv), _buf_rv);

							if (_buf.empty() || _buf_rv.empty()){
								std::cout << "wrong\n";
								std::cout << _buf << "\n";
								std::cout << _buf_rv << "\n\n";

								std::cout << "create: err reading paired fasta files:\n"
										<< filelist_[i].filename << "\n"
										<< filelist_[i + 1].filename
										<< " ... exit \n";
								exit(1);
							} // if

						}
					} // if
				} // while

				filelist_[i].numReads = numReads;
				filelist_[i + 1].numReads = numReads;

				++i;

			} else { // unpaired file

				while (!std::getline((*fH), _buf).eof()) {
					++ numReads;
					if (!_buf.empty()) {
						if (_buf.at(0) == '@') { //fastq
							_is_pos = fH->tellg();
							_is_pos -= (_buf.size() + 1);
							_rIdx.push_back(ReadIdx(i, _is_pos));
							for (int j = 0; j < 3; ++j) {
								std::getline((*fH), _buf);
								if (_buf.empty()) {
									std::cout << "create: err reading fastq file:\n"
									<< filelist_[i].filename
									<< " ... exit \n";
									exit(1);
								}

							}
						} else if (_buf.at(0) == '>') { //fasta
							_is_pos = fH->tellg();
							_is_pos -= (_buf.size() + 1);
							_rIdx.push_back(ReadIdx(i, _is_pos));
							std::getline((*fH), _buf);
							if (_buf.empty()){
								std::cout << "create: err reading fasta file:\n"
									  << filelist_[i].filename
									  << " ... exit \n";
								exit(1);
							} // if
						} // else if
					} // if
				} // while

				filelist_[i].numReads = numReads;

			} // else
		} // for
	} // create


};

#endif /* READINDEXER_H_ */
