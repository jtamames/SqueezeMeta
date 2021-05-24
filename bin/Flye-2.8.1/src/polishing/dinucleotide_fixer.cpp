//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)
//
#include "dinucleotide_fixer.h"
#include "alignment.h"

void DinucleotideFixer::fixBubble(Bubble& bubble) const
{
	auto likelihood = [this](const std::string& candidate, 
						     const std::vector<std::string>& branches)
	{
		Alignment align(branches.size(), _subsMatrix);
		AlnScoreType score = align.globalAlignment(candidate, branches);
		return score;
	};

	auto runPair = this->getDinucleotideRuns(bubble.candidate);
	if (runPair.second < 3) return;
	
	//try to increase / decrease dinucleotide polimer len and see what happens
	std::string increased = bubble.candidate;
	increased.insert(runPair.first, bubble.candidate.substr(runPair.first, 2));

	std::string decreased = bubble.candidate;
	decreased.erase(runPair.first, 2);

	AlnScoreType normalScore = likelihood(bubble.candidate, bubble.branches);
	AlnScoreType increasedScore = likelihood(increased, bubble.branches);
	AlnScoreType decreasedScore = likelihood(decreased, bubble.branches);

	/*
	if (increasedScore > normalScore || decreasedScore > normalScore)
	{
		std::cerr << bubble.candidate << std::endl << increased 
				  << std::endl << decreased << std::endl;
		std::cerr << normalScore << " " 
				  << increasedScore << " " << decreasedScore << std::endl << std::endl;
	}*/

	if (increasedScore > normalScore)
	{
		bubble.candidate = increased;
		StepInfo info;
		info.sequence = increased;
		bubble.polishSteps.push_back(info);
	}
	else if (decreasedScore > normalScore)
	{
		bubble.candidate = decreased;
		StepInfo info;
		info.sequence = decreased;
		bubble.polishSteps.push_back(info);
	}
}

std::pair<int, int>
DinucleotideFixer::getDinucleotideRuns(const std::string& sequence) const
{
	int maxRun = 0;
	int maxPos = 0;
	if (sequence.length() < 2) return {0, 0};

	for (size_t shift = 0; shift < 2; ++shift)
	{
		std::string prevDinuc = "--";
		int curRun = 0;
		for (size_t pos = 0; pos < (sequence.length() - shift) / 2; ++pos)
		{
			std::string curDinuc = sequence.substr(pos * 2 + shift, 2);
			//std::cerr << curDinuc << " ";
			if (curDinuc != prevDinuc)
			{
				if (curRun > maxRun && prevDinuc[0] != prevDinuc[1])
				{
					maxRun = curRun;
					maxPos = pos * 2 + shift - curRun * 2;
				}
				curRun = 1;
				prevDinuc = curDinuc;
			}
			else
			{
				++curRun;
			}
		}
		//std::cerr << std::endl;

	}
	
	/*if (maxRun > 2)
	{
		std::cerr << sequence << std::endl;
		std::cerr << maxPos << " " << maxRun << std::endl;
		std::cerr << std::endl;
	}*/

	return {maxPos, maxRun};
}
