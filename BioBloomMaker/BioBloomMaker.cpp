/*
 * BioBloom.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include "BloomFilterGenerator.h"
#include "BloomFilter/BloomFilterInfo.hpp"
#include "Common/SeqEval.h"
#include <boost/unordered/unordered_map.hpp>
#include <getopt.h>
#include "config.h"
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "biobloommaker"

namespace opt {
/** The number of parallel threads. */
static unsigned threads = 1;
}

void printVersion() {
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
	"Written by Justin Chu.\n"
	"\n"
	"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	static const char dialog[] =
		"Usage: biobloommaker -p [FILTERID] [OPTION]... [FILE]...\n"
		"Usage: biobloommaker -p [FILTERID] -r 0.2 [FILE]... [FASTQ1] [FASTQ2] \n"
		"Creates a bf and txt file from a list of fasta files. The input sequences are\n"
		"cut into a k-mers with a sliding window and their hash signatures are inserted\n"
		"into a bloom filter.\n"
		"\n"
		"  -p, --file_prefix=N    Filter prefix and filter ID. Required option.\n"
		"  -o, --output_dir=N     Output location of the filter and filter info files.\n"
		"  -h, --help             Display this dialog.\n"
		"  -v  --version          Display version information.\n"
		"  -t, --threads=N        The number of threads to use. Experimental. [1]\n"
		"  -i, --inclusive        If one paired read matches, both reads will be included\n"
		"                         in the filter. \n"
		"\nAdvanced options:\n"
		"  -f, --fal_pos_rate=N   Maximum false positive rate to use in filter. [0.0075]\n"
		"  -g, --hash_num=N       Set number of hash functions to use in filter instead\n"
		"                         of automatically using calculated optimal number of\n"
		"                         functions.\n"
		"  -k, --kmer_size=N      K-mer size to use to create filter. [25]\n"
		"  -s, --subtract=N       Path to filter that you want to uses to prevent the\n"
		"                         addition of k-mers contained into new filter. You may\n"
		"                         only use filters with k-mer sizes equal the one you\n"
		"                         wish to create.\n"
		"  -n, --num_ele=N        Set the number of expected elements. If set to 0 number\n"
		"                         is determined from sequences sizes within files. [0]\n"
		"  -r, --progressive=N    Progressive filter creation. The score threshold is\n"
		"                         specified by N, which may be either a floating point score\n"
		"                         between 0 and 1 or a positive integer.  If N is a\n"
		"                         positive integer, it is interpreted as the minimum\n"
		"                         number of contiguous matching bases required for a\n"
		"                         match.\n"
		"\n"
		"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

int main(int argc, char *argv[]) {

	bool die = false;

	//switch statement variable
	int c;

	//command line variables
	double fpr = 0.0075;
	string filterPrefix = "";
	string outputDir = "";
	unsigned kmerSize = 25;
	unsigned hashNum = 0;
	string subtractFilter = "";
	size_t entryNum = 0;
	double progressive = -1;
	bool inclusive = false;
	SeqEval::EvalMode evalMode = SeqEval::EVAL_STANDARD;

	//long form arguments
	static struct option long_options[] = {
			{
					"fal_pos_rate", required_argument, NULL, 'f' }, {
					"file_prefix", required_argument, NULL, 'p' }, {
					"output_dir", required_argument, NULL, 'o' }, {
					"threads", required_argument, NULL, 't' }, {
					"inclusive", no_argument, NULL, 'i' }, {
					"version", no_argument, NULL, 'v' }, {
					"hash_num", required_argument, NULL, 'g' }, {
					"kmer_size", required_argument, NULL, 'k' }, {
					"subtract",	required_argument, NULL, 's' }, {
					"num_ele", required_argument, NULL, 'n' }, {
					"help", no_argument, NULL, 'h' }, {
					"progressive", required_argument, NULL, 'r' }, {
					NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:p:o:k:n:g:hvs:n:t:r:i", long_options,
			&option_index)) != -1) {
		switch (c) {
		case 'f': {
			stringstream convert(optarg);
			if (!(convert >> fpr)) {
				cerr << "Error - Invalid set of bloom filter parameters! f: "
						<< optarg << endl;
				return 0;
			}
			if (fpr > 1) {
				cerr << "Error -f cannot be greater than 1 " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			filterPrefix = optarg;
			break;
		}
		case 'o': {
			outputDir = optarg;
			if (outputDir.at(outputDir.length() - 1) != '/') {
				outputDir = outputDir + '/';
			}
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter! t: " << optarg << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'i': {
			inclusive = true;
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> kmerSize)) {
				cerr << "Error - Invalid set of bloom filter parameters! k: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> hashNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'n': {
			stringstream convert(optarg);
			if (!(convert >> entryNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! n: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'v': {
			printVersion();
			break;
		}
		case 's': {
			stringstream convert(optarg);
			if (!(convert >> subtractFilter)) {
				cerr << "Error - Invalid set of bloom filter parameters! s: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			unsigned matchLen;
			// if arg is a positive integer > 1, interpret as minimum match
			// length in bases
			if ((convert >> matchLen) && matchLen > 1) {
				progressive = matchLen;
				evalMode = SeqEval::EVAL_MIN_MATCH_LEN;
			} else {
				// not a positive integer > 1, so interpret as floating
				// point score between 0 and 1
				stringstream convert2(optarg);
				if (!(convert2 >> progressive)) {
					cerr << "Error - Invalid set of bloom filter parameters! r: "
						<< optarg << endl;
					return 0;
				}
				if (progressive < 0 || progressive > 1) {
					cerr << "Error - s must be a positive integer or a floating "
						<< "point between 0 and 1. Input given:"
						<< optarg << endl;
					exit(EXIT_FAILURE);
				}
				evalMode = SeqEval::EVAL_STANDARD;
			}
			break;
		}
		default: {
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
	omp_set_num_threads(opt::threads);
#endif

	//Stores fasta input file names
	vector<string> inputFiles;

	while (optind < argc) {
		inputFiles.push_back(argv[optind]);
		optind++;
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Need Input File" << endl;
		die = true;
	}
	if (filterPrefix.size() == 0) {
		cerr << "Need Filter Prefix ID" << endl;
		die = true;
	}
	if (filterPrefix.find('/') != string::npos) {
		cerr << "Prefix ID cannot have '/' characters" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//set number of hash functions used
	if (hashNum == 0) {
		//get optimal number of hash functions
		hashNum = unsigned(-log(fpr) / log(2));
	}

	string file1 = "";
	string file2 = "";

	if (progressive != -1) {
		if (inputFiles.size() > 2) {
			file1 = inputFiles.back();
			inputFiles.pop_back();
			file2 = inputFiles.back();
			inputFiles.pop_back();
			cerr << "Building Bloom filter in progessive mode. ";
			switch(evalMode) {
				case SeqEval::EVAL_MIN_MATCH_LEN:
					cerr << "Min match length = "
						<< (unsigned)round(progressive)
						<< " bp" << endl;
					break;
				case SeqEval::EVAL_STANDARD:
				default:
					cerr << "Score threshold = "
						<< progressive << endl;
					break;
			}
		} else {
			cerr << "require a least 3 input when using progressive mode"
					<< endl;
			exit(1);
		}
	}

	//create filter
	BloomFilterGenerator filterGen(inputFiles, kmerSize, hashNum, entryNum);

	if (entryNum == 0) {
		filterGen = BloomFilterGenerator(inputFiles, kmerSize, hashNum);
		entryNum = filterGen.getExpectedEntries();
	}

	BloomFilterInfo info(filterPrefix, kmerSize, hashNum, fpr, entryNum,
			inputFiles);

	//get calculated size of Filter
	size_t filterSize = info.getCalcuatedFilterSize();
	cerr << "Allocating " << filterSize << " bits of space for filter and will output filter this size" << endl;
	filterGen.setFilterSize(filterSize);

	size_t redundNum = 0;
	//output filter
	if (!subtractFilter.empty() && progressive != -1) {
		createMode mode = PROG_STD;
		if (inclusive) {
			mode = PROG_INC;
		}
		redundNum = filterGen.generateProgressive(
				outputDir + filterPrefix + ".bf", progressive, file1, file2,
				mode, evalMode, subtractFilter);
	} else if (!subtractFilter.empty()) {
		redundNum = filterGen.generate(outputDir + filterPrefix + ".bf",
				subtractFilter);
	} else if (progressive != -1) {
		createMode mode = PROG_STD;
		if (inclusive) {
			mode = PROG_INC;
		}
		redundNum = filterGen.generateProgressive(
				outputDir + filterPrefix + ".bf", progressive, file1, file2,
				mode, evalMode);
	} else {
		redundNum = filterGen.generate(outputDir + filterPrefix + ".bf");
	}
	info.setTotalNum(filterGen.getTotalEntries());
	info.setRedundancy(redundNum);

	//code for redundancy checking
	//calculate redundancy rate
	double redunRate = double(redundNum) / double(filterGen.getTotalEntries())
			- info.getRedundancyFPR();
	if (redunRate > 0.25) {
		cerr
				<< "The ratio between redundant k-mers and unique k-mers is approximately: "
				<< redunRate << endl;
		cerr
				<< "Consider checking your files for duplicate sequences and adjusting them accordingly.\n"
						"High redundancy will cause filter sizes used overestimated, potentially resulting in a larger than needed filter.\n"
						"Alternatively you can set the number of elements wanted in the filter with (-n) and ignore this message."
				<< endl;
	}

	//output info
	info.printInfoFile(outputDir + filterPrefix + ".txt");
	cerr << "Filter Creation Complete." << endl;

	return 0;
}
