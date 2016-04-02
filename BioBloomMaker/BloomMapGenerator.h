/*
 * BloomMapGenerator.h
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#ifndef BLOOMMAPGENERATOR_H_
#define BLOOMMAPGENERATOR_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include <google/dense_hash_map>
#include "bloomfilter/BloomMapSS.hpp"
#include "bloomfilter/RollingHashIterator.h"
#include "Common/Options.h"

using namespace std;

class BloomMapGenerator {
public:
	explicit BloomMapGenerator(vector<string> const &filenames,
			unsigned kmerSize, size_t numElements);

	void generate(const string &filePrefix, double fpr);

	virtual ~BloomMapGenerator();
private:

	unsigned m_kmerSize;
//	unsigned m_hashNum;
	size_t m_expectedEntries;
	size_t m_totalEntries;
	vector<string> m_fileNames;
	google::dense_hash_map<ID,string> m_headerIDs;

	//helper methods
	//TODO: collision detection
	inline size_t loadSeq(BloomMapSS<ID> &bloomMap, const string& seq, ID value) {
		if (seq.size() < m_kmerSize)
			return 0;
		size_t count = 0;
		/* init rolling hash state and compute hash values for first k-mer */
		RollingHashIterator itr(seq, m_kmerSize, bloomMap.getSeedValues());
		while (itr != itr.end()) {
			vector<ID> temp = bloomMap.at(*itr);
//			bool collision = false;
//			for (unsigned i = 0; i < temp.size(); ++i) {
//				if (temp[i] != 0 && temp[i] != value) {
//					collision = true;
//				}
//			}
//			if (collision) {
//				cerr << value;
//				for (unsigned i = 0; i < temp.size(); ++i) {
//					if (temp[i] != 0) {
//						cerr << "\t" << temp[i];
//					} else {
//						cerr << "\t-";
//					}
//				}
//				cerr << endl;
//			}
			count += !bloomMap.insertAndCheck(*itr, value);
			++itr;
		}
		return count;
	}

	inline void writeIDs(const string& filename, google::dense_hash_map<ID,string> headerIDs) {
		std::ofstream file;
		file.open(filename.c_str());
		assert(file);
		for (google::dense_hash_map<ID,string>::iterator itr = headerIDs.begin(); itr != headerIDs.end(); ++itr) {
			file << (*itr).first << "\t" << (*itr).second << "\n";
			assert(file);
		}
		file.close();
	}

		
};

#endif /* BLOOMMAPGENERATOR_H_ */
