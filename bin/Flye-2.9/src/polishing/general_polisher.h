//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "bubble.h"
#include "subs_matrix.h"
#include "alignment.h"

class GeneralPolisher
{
public:
	GeneralPolisher(const SubstitutionMatrix& subsMatrix):
		_subsMatrix(subsMatrix)
	{}
	void polishBubble(Bubble& bubble) const;

private:
	StepInfo makeStep(const std::string& candidate, 
					  const std::vector<std::string>& branches,
					  Alignment& align) const;

	const SubstitutionMatrix& _subsMatrix;
};
