cmake_minimum_required(VERSION 3.8)
project(biobloom)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        BioBloomCategorizer/BioBloomCategorizer.cpp
        BioBloomCategorizer/BioBloomClassifier.cpp
        BioBloomCategorizer/BioBloomClassifier.h
        BioBloomCategorizer/MultiFilter.cpp
        BioBloomCategorizer/MultiFilter.h
        BioBloomCategorizer/ResultsManager.cpp
        BioBloomCategorizer/ResultsManager.h
        BioBloomMaker/BioBloomMaker.cpp
        BioBloomMaker/BloomFilterGenerator.cpp
        BioBloomMaker/BloomFilterGenerator.h
        BioBloomMaker/WindowedFileParser.cpp
        BioBloomMaker/WindowedFileParser.h
        Common/BloomFilter.cpp
        Common/BloomFilter.h
        Common/BloomFilterInfo.cpp
        Common/BloomFilterInfo.h
        Common/city.cc
        Common/city.h
        Common/citycrc.h
        Common/Dynamicofstream.cpp
        Common/Dynamicofstream.h
        Common/Fcontrol.cpp
        Common/Fcontrol.h
        Common/gzstream.C
        Common/gzstream.h
        Common/IOUtil.h
        Common/Options.cpp
        Common/Options.h
        Common/ReadsProcessor.cpp
        Common/ReadsProcessor.h
        Common/SeqEval.h
        Common/Sequence.cpp
        Common/Sequence.h
        Common/SignalHandler.cpp
        Common/SignalHandler.h
        Common/StringUtil.h
        Common/Uncompress.cpp
        Common/Uncompress.h
        DataLayer/FastaIndex.h
        DataLayer/FastaInterleave.h
        DataLayer/FastaReader.cpp
        DataLayer/FastaReader.h
        DataLayer/FastaWriter.cpp
        DataLayer/FastaWriter.h
        DataLayer/kseq.h
        DataLayer/Options.h
        Tests/BloomFilterCategorizerTests.cpp
        Tests/BloomFilterInfoTests.cpp
        Tests/BloomFilterMakerTests.cpp
        Tests/BloomFilterTests.cpp
        Tests/ReadProcessorTests.cpp
        Tests/SeqEvalTests.cpp
        Tests/WindowedParserTests.cpp)

add_executable(biobloom ${SOURCE_FILES})