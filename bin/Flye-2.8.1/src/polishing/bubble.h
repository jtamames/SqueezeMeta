//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>

#include "subs_matrix.h"

struct StepInfo 
{
	std::string sequence;
	AlnScoreType score;

	StepInfo(): score(0.0f) {}
};

struct Bubble
{
	std::string header;
	int position;

	std::string candidate;
	std::vector<std::string> branches;
	std::vector<StepInfo> polishSteps;
};
