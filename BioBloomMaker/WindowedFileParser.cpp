/*
 * WindowedFileParser.cpp
 *
 * Currently for fasta files only and needs fasta index. Parse file line by line.
 *
 *  Created on: Jul 18, 2012
 *      Author: cjustin
 */
#include "WindowedFileParser.h"
#include <sstream>
#include "DataLayer/FastaIndex.h"

WindowedFileParser::WindowedFileParser(string const &fileName,
		unsigned windowSize):
		m_windowSize(windowSize)
{
	m_fastaFileHandle.open(fileName.c_str(), ifstream::in);
	assert(m_fastaFileHandle);

	//create in memory index
	WindowedFileParser::initializeIndex(fileName);
	m_currentHeader = "";
	setLocationByHeader(m_headers[0]);
}

const vector<string> WindowedFileParser::getHeaders() const
{
	return m_headers;
}

//sets the location in the file to the start of the sequence given a header
void WindowedFileParser::setLocationByHeader(string const &header)
{
	m_sequenceNotEnd = true;
	m_currentHeader = header;
	m_fastaFileHandle.seekg(m_fastaIndex[header].start, ios::beg);
	m_currentCharNumber = 0;
	string bufferString;
	getline(m_fastaFileHandle, m_currentString);
	while ((m_currentString.length() < m_windowSize)
			&& (m_currentCharNumber < m_fastaIndex[m_currentHeader].size)
			&& getline(m_fastaFileHandle, bufferString))
	{
		m_currentString += bufferString;
	}
	m_nextNonATCG = ReadsProcessor::prepSeq(m_currentString);
	m_reset = true;
	m_nextNonATCG = extendUntilNonATCG(m_nextNonATCG);
	m_currentLinePos = 0;
}

size_t WindowedFileParser::getSequenceSize(string const &header) const
{
	return m_fastaIndex.at(header).size;
}

/*
 * Extends until a long enough stretch of Non ATCG is encountered
 */
unsigned WindowedFileParser::extendUntilNonATCG(unsigned nextNonATCG){
	//base case
	if (nextNonATCG >= m_windowSize + m_currentLinePos) {
		return nextNonATCG;
	}

	m_reset = true;
	m_currentString.erase(0, nextNonATCG + 1);
	m_currentLinePos = 0;
	while ((m_currentString.length() < m_windowSize)
			&& (m_currentCharNumber < m_fastaIndex[m_currentHeader].size)
			&& getline(m_fastaFileHandle, m_bufferString)) {
		m_currentString += m_bufferString;
	}
	nextNonATCG = ReadsProcessor::prepSeq(m_currentString);

	if (m_currentString.length() < m_windowSize) {
		m_sequenceNotEnd = false;
	}
	return(extendUntilNonATCG(nextNonATCG));
}


/*
 * Return the next string in sliding window, also cleans and formats
 * sequences using ReadProcessor
 */
char* WindowedFileParser::getNextSeq()
{
	char* outputStr = m_currentString.substr(m_currentLinePos++, m_windowSize).c_str();
	m_nextNonATCG = extendUntilNonATCG(m_nextNonATCG);
	if (m_currentString.length() < m_windowSize + m_currentLinePos) {
		m_currentString.erase(0, m_currentLinePos);
		m_currentLinePos = 0;
		//grow the sequence to match the correct window size
		//stop if there are no more lines left in fasta file
		while (m_fastaFileHandle.is_open()
				&& (m_currentString.length() < m_windowSize)
				&& (m_currentCharNumber < m_fastaIndex[m_currentHeader].size)
				&& getline(m_fastaFileHandle, m_bufferString)) {
			m_currentString += m_bufferString;
		}
		m_nextNonATCG = ReadsProcessor::prepSeq(m_currentString);

		//if there is not enough sequence for a full kmer
		if (m_currentString.length() < m_windowSize) {
			m_sequenceNotEnd = false;
		}
	}
	return outputStr;
}

/*
 * Obtains out and in char for rolling hash
 * returns if a reset was performed.
 */
bool WindowedFileParser::getNextChar(char &out, char &in){
	m_reset = false;
	m_nextNonATCG = extendUntilNonATCG(m_nextNonATCG);
	if(m_reset){
		return m_reset;
	}
	out = m_currentString[m_currentLinePos];
	in = m_currentString[m_currentLinePos++ + m_windowSize];
	if (m_nextNonATCG < m_windowSize + m_currentLinePos) {
		m_currentString.erase(0, m_currentLinePos);
		m_currentLinePos = 0;
		//grow the sequence to match the correct window size
		//stop if there are no more lines left in fasta file
		while (m_fastaFileHandle.is_open()
				&& (m_currentString.length() < m_windowSize)
				&& (m_currentCharNumber < m_fastaIndex[m_currentHeader].size)
				&& getline(m_fastaFileHandle, m_bufferString))
		{
			m_currentString += m_bufferString;
		}
		m_nextNonATCG = ReadsProcessor::prepSeq(m_currentString);

		//if there is not enough sequence for a full kmer
		if (m_currentString.length() < m_windowSize) {
			m_sequenceNotEnd = false;
		}
	}
	return m_reset;
}

bool WindowedFileParser::notEndOfSeqeunce() const
{
	return m_sequenceNotEnd;
}

/*
 * Initializes fasta index in memory.
 * If no index exists derive one from input file
 * Input file refers to input fasta file not the index file.
 */
void WindowedFileParser::initializeIndex(string const &fileName)
{
	string faiFile = fileName + ".fai";
	//look for fasta index file
	//append .fai to end of filename
	ifstream indexFile;
	indexFile.open(faiFile.c_str(), ifstream::in);
	if (!indexFile) {
		cerr << "Fasta files must be indexed. Use samtools faidx." << endl;
		exit(1);
	}

	string line;
	//read file and populate index struct
	if (indexFile.is_open()) {
		while (getline(indexFile, line)) {
			stringstream ss(line);
			string header;
			FastaIndexValue value;
			ss >> header >> value.size >> value.start >> value.bpPerLine
					>> value.charsPerLine;
			value.index = m_headers.size();
			m_fastaIndex[header] = value;
			m_headers.push_back(header);
		}
		indexFile.close();
	}
}

WindowedFileParser::~WindowedFileParser()
{
	m_fastaFileHandle.close();
}

