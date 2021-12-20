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

	//Logger::get().debug() << "Stitching " << path.name;

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
			//Logger::get().debug() << "\tPiece " << sequence.length() << " " 
			//	<< leftCut << " " << rightCut << " " << rightCut - path.overlaps[i].curBegin;
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
	//typedef std::pair<const ContigPath*, size_t> AlnTask;
	struct AlnTask
	{
		const ContigPath* path;
		size_t contigId;
		OverlapRange adjustedOverlap;
	};

	AlignmentsMap alnMap;
	std::mutex mapMutex;
	std::function<void(const AlnTask&)> alnFunc =
	[&alnMap, &mapMutex](const AlnTask& task)
	{
		const ContigPath* path = task.path;
		size_t i = task.contigId;
		auto curOverlap = task.adjustedOverlap;

		std::vector<CigOp> cigar;
		const float maxErr = 0.3;
		std::string alignedLeft;
		std::string alignedRight;
		getAlignmentCigarKsw(path->sequences[i], curOverlap.curBegin, curOverlap.curRange(),
			   			     path->sequences[i + 1], curOverlap.extBegin, curOverlap.extRange(),
			   			   	 maxErr, cigar);
		decodeCigar(cigar, path->sequences[i], curOverlap.curBegin,
				 	path->sequences[i + 1], curOverlap.extBegin,
				 	alignedLeft, alignedRight);

		{
			std::lock_guard<std::mutex> lock(mapMutex);
			alnMap[&path->overlaps[i]] = {alignedLeft, alignedRight, 
						 				  curOverlap.curBegin, 
										  curOverlap.extBegin};
		}
	};

	std::vector<AlnTask> tasks;
	for (auto& path : contigs)
	{
		int32_t prevSwitch = 0;
		for (size_t i = 0; i < path.sequences.size() - 1; ++i)
		{
			OverlapRange curOverlap = path.overlaps[i];

			//don't compute alignment for regions we know will
			//not be used for stitching
			int32_t beginShift = prevSwitch - curOverlap.curBegin;
			if (beginShift > 0 && 
				beginShift < std::min(curOverlap.curRange(), curOverlap.extRange()))
			{
				curOverlap.curBegin += beginShift;
				curOverlap.extBegin += beginShift;
			}
			prevSwitch = curOverlap.extBegin;

			//in case of long reads, only consider last 20k of the overlap to
			//save memory during pairwise alignmemnt
			const int32_t MAX_ALIGNMENT = 20000;
			int32_t endShift = std::min(curOverlap.curRange(), 
										curOverlap.extRange()) - MAX_ALIGNMENT;
			if (endShift > 0)
			{
				curOverlap.curEnd -= endShift;
				curOverlap.extEnd -= endShift;
			}

			tasks.push_back({&path, i, curOverlap});
		}
	}
	processInParallel(tasks, alnFunc, Parameters::get().numThreads, verbose);

	return alnMap;
}


std::pair<int32_t, int32_t> 
ConsensusGenerator::getSwitchPositions(const AlignmentInfo& aln,
									   int32_t prevSwitch)
{
	const int MIN_SEGMENT = 500;
	const int MIN_MATCH = 15;

	int leftPos = aln.startOne;
	int rightPos = aln.startTwo;
	int matchRun = 0;
	for (size_t i = 0; i < aln.alnOne.length(); ++i)
	{
		if (aln.alnOne[i] != '-') ++leftPos;
		if (aln.alnTwo[i] != '-') ++rightPos;

		if (aln.alnOne[i] != '-' && aln.alnTwo[i] != '-' &&
			leftPos > prevSwitch + MIN_SEGMENT)
		{
			++matchRun;
		}
		else
		{
			matchRun = 0;
		}
		if (matchRun == MIN_MATCH)
		{
			return {leftPos, rightPos};
		}
	}

	//Logger::get().debug() << "Prev: " << prevSwitch;
	//for (size_t i = 0; i < aln.alnOne.size() / 100; ++i)
	//{
	//	Logger::get().debug() << aln.alnOne.substr(i * 100, 100);
	//	Logger::get().debug() << aln.alnTwo.substr(i * 100, 100);
	//	Logger::get().debug() << "";
	//}

	//Logger::get().info() << "No jump found!";
	prevSwitch = std::max(prevSwitch + 1, aln.startOne);
	return {prevSwitch, aln.startTwo};
}
