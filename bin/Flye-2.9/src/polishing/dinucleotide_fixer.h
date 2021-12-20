//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "subs_matrix.h"
#include "bubble.h"

class DinucleotideFixer
{
public:
	DinucleotideFixer(const SubstitutionMatrix& subsMatrix):
		_subsMatrix(subsMatrix)
	{}
	void fixBubble(Bubble& bubble) const;

private:
	std::pair<int, int> getDinucleotideRuns(const std::string& sequence) const;

	const SubstitutionMatrix& _subsMatrix;
};
