//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "consensus_generator.h"
#include "../sequence/alignment.h"

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/matrix.h"


std::vector<FastaRecord> 
	ConsensusGenerator::generateConsensuses(const std::vector<ContigPath>& contigs, 
											bool verbose)
{
	if (verbose) Logger::get().info() << "Generating sequence";
	std::vector<std::vector<AlignmentInfo>> allAlignments;
	std::vector<FastaRecord> consensuses;

	auto alnMap = this->generateAlignments(contigs, verbose);
	//then, generate contig sequences
	for (size_t i = 0; i < contigs.size(); ++i)
	{
		if (contigs[i].sequences.empty()) continue;
		if (contigs[i].sequences.size() == 1)
		{
			FastaRecord rec(contigs[i].sequences.front(), contigs[i].name, 
							FastaRecord::ID_NONE);
			consensuses.push_back(rec);
		}
		else
		{
			consensuses.push_back(this->generateLinear(contigs[i], alnMap));
		}
	}
	return consensuses;
}


FastaRecord ConsensusGenerator::generateLinear(const ContigPath& path, 
											   const AlignmentsMap& alnMap)
{
	std::vector<FastaRecord> contigParts;

	auto prevSwitch = std::make_pair(0, 0);
	std::string contigSequence;
	for (size_t i = 0; i < path.sequences.size(); ++i)
	{
		auto& sequence = path.sequences[i];
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = sequence.length();
		if (i != path.sequences.size() - 1)
		{
			auto curSwitch = 
				this->getSwitchPositions(alnMap.at(&path.overlaps[i]),
										 prevSwitch.second);
			rightCut = curSwitch.first;
			prevSwitch = curSwitch;
		}

		if (rightCut - leftCut > 0)	//shoudn't happen, but just in case
		{
			contigSequence += sequence.substr(leftCut, rightCut - leftCut).str();
		}
	}
	int32_t cutLen = contigSequence.length() - (path.trimLeft + path.trimRight);
	if (cutLen > 0)
	{
		contigSequence = contigSequence.substr(path.trimLeft, cutLen);
	}
	return FastaRecord(DnaSequence(contigSequence), path.name, 
					   FastaRecord::ID_NONE);
}


ConsensusGenerator::AlignmentsMap 
	ConsensusGenerator::generateAlignments(const std::vector<ContigPath>& contigs,
										   bool verbose)
{
	typedef std::pair<const ContigPath*, size_t> AlnTask;

	AlignmentsMap alnMap;
	std::mutex mapMutex;
	std::function<void(const AlnTask&)> alnFunc =
	[&alnMap, &mapMutex](const AlnTask& task)
	{
		const ContigPath* path = task.first;
		size_t i = task.second;

		std::vector<CigOp> cigar;
		const float maxErr = 0.3;
		std::string alignedLeft;
		std::string alignedRight;
		getAlignmentCigarKsw(path->sequences[i], path->overlaps[i].curBegin, path->overlaps[i].curRange(),
			   			     path->sequences[i + 1], path->overlaps[i].extBegin, path->overlaps[i].extRange(),
			   			   	 maxErr, cigar);
		decodeCigar(cigar, path->sequences[i], path->overlaps[i].curBegin,
				 	path->sequences[i + 1], path->overlaps[i].extBegin,
				 	alignedLeft, alignedRight);

		{
			std::lock_guard<std::mutex> lock(mapMutex);
			alnMap[&path->overlaps[i]] = {alignedLeft, alignedRight, 
						 				  path->overlaps[i].curBegin, 
										  path->overlaps[i].extBegin};
		}
	};

	std::vector<AlnTask> tasks;
	for (auto& path : contigs)
	{
		for (size_t i = 0; i < path.sequences.size() - 1; ++i)
		{
			tasks.emplace_back(&path, i);
		}
	}
	processInParallel(tasks, alnFunc, Parameters::get().numThreads, verbose);

	return alnMap;
}


std::pair<int32_t, int32_t> 
ConsensusGenerator::getSwitchPositions(const AlignmentInfo& aln,
									int32_t prevSwitch)
{
	int leftPos = aln.startOne;
	int rightPos = aln.startTwo;
	int matchRun = 0;
	for (size_t i = 0; i < aln.alnOne.length(); ++i)
	{
		if (aln.alnOne[i] != '-') ++leftPos;
		if (aln.alnTwo[i] != '-') ++rightPos;

		if (aln.alnOne[i] == aln.alnTwo[i] &&
			leftPos > prevSwitch + Config::get("maximum_jump"))
		{
			++matchRun;
		}
		else
		{
			matchRun = 0;
		}
		if (matchRun == (int)Parameters::get().kmerSize)
		{
			return {leftPos, rightPos};
		}
	}

	//Logger::get().info() << "No jump found!";
	prevSwitch = std::max(prevSwitch + 1, aln.startOne);
	return {prevSwitch, aln.startTwo};
}
