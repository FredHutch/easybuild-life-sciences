//============================================================================
// Project     : Diversifier
// Name        : file_manip.hpp
// Author      : Xiao Yang
// Created on  : Aug 22, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef FILE_MANIP_HPP_
#define FILE_MANIP_HPP_

#include <iostream>
#include <string>
#include <dirent.h>
#include <fstream>

namespace xny{

	template <typename Output> bool dir_list(const std::string& name, Output out) {
		DIR* dir = opendir(name.c_str());
		if (dir == 0) return false;

		dirent* entry = readdir(dir);
		while (entry != 0) {
		  *(out++) = entry->d_name;
		  entry = readdir(dir);
		}
		closedir(dir);
	  return true;
	} // dir_list

	/*
	 * find the last/first [delim] in [str] if [dir] is set to be true/false
	 */
	inline std::string get_suffix (const std::string& str, bool last,
			const std::string& delim) {
		std::string suf = str;
		int pos = suf.find_last_of(delim);
		if (!last) pos = suf.find_first_of (delim);
		if (pos != (int) std::string::npos)
			return suf.substr(pos + 1, (int) suf.length() - pos);
		else return "";
	}

	/*
	 * find the first/last [delim] in [str] if [dir] is set to be true/false
	 */
	inline std::string get_prefix (const std::string& str, bool first,
			const std::string& delim) {
		std::string prfx = str;
		int pos = prfx.find_last_of(delim);
		if (first) pos = prfx.find_first_of(delim);
		if (pos!= (int) std::string::npos)
			return prfx.substr(0, pos);
		else return "";
	}

	inline void openfile (std::ofstream& outfile, const std::string& fileName){
		if (!outfile.is_open()) {
			outfile.open(fileName.c_str());
			//outfile.open(fileName.c_str(), std::fstream::app);
			if (! outfile.good() ){
				std::cout << "[WARNING] opening file: " << fileName << " failed!\n"
						  << "Will write to standard output ... \n";
			}
		}
	}

	inline void closefile(std::ofstream& outfile){
		if (outfile.good())	outfile.close();
	}

/*
	inline std::string get_suf_filename (const std::string& filePath) {
		std::string suf = filePath;
		int dotPos = suf.find_last_of(".");
		if (dotPos != (int) std::string::npos)
			return suf.substr(dotPos + 1, (int) suf.length() - dotPos);
		return std::string ("");
	}

	inline std::string get_pref_filename (const std::string& filePath) {
		std::string prfx = filePath;
		int dotPos = prfx.find_last_of(".");
		if (dotPos!= (int) std::string::npos) {
			int slashPos = prfx.find_last_of("/");
			if (slashPos != (int) std::string::npos){
				prfx = prfx.substr(slashPos, dotPos-slashPos);
			}
			else prfx = prfx.substr(0, dotPos);
		}
		return prfx;
	}

	inline std::string get_filename_frompath (const std::string& filePath){
		std::string suf = filePath;
		int slashPos = filePath.find_last_of("/");
		if (slashPos != (int) std::string::npos){
			return suf.substr(slashPos + 1, (int) suf.length() - slashPos);
		}
		return filePath;
	}
*/
}

#endif /* FILE_MANIP_HPP_ */
