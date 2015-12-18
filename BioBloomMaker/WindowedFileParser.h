/*
 * WindowedFileParser.h
 *
 *  Created on: Jul 18, 2012
 *      Author: cjustin
 */

#ifndef WINDOWEDFILEPARSER_H_
#define WINDOWEDFILEPARSER_H_
#include <vector>
#include <boost/unordered/unordered_map.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include "DataLayer/FastaReader.h"
#include <deque>

using namespace std;
using namespace boost;

class WindowedFileParser {
public:
	//constructor/destructor
	explicit WindowedFileParser(const string &fileName, unsigned windowSize);
	const vector<string> getHeaders() const;
	void setLocationByHeader( const string &header);
	size_t getSequenceSize( const string &header) const;
	char* getNextSeq();
	bool getNextChar(char &out, char &in);

	bool notEndOfSeqeunce() const;

	virtual ~WindowedFileParser();

private:
	struct FastaIndexValue {
		size_t index;
		size_t size;
		size_t start;
		size_t bpPerLine;
		size_t charsPerLine;
	};

	unordered_map<string, FastaIndexValue> m_fastaIndex;
	ifstream m_fastaFileHandle;
	unsigned m_windowSize;
	vector<string> m_headers;
	string m_currentHeader;
	unsigned m_currentCharNumber;
	unsigned m_currentLinePos;
	string m_window;
	string m_currentString;
	bool m_sequenceNotEnd;

	unsigned m_nextNonATCG;
	bool m_reset;

	string m_bufferString; //so reallocation does not need to occur

	//helper methods
	void initializeIndex(string const &fileName);

	unsigned extendUntilNonATCG(unsigned nextNonATCG);
};

#endif /* WINDOWEDFILEPARSER_H_ */
