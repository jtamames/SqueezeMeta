//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include "../common/matrix.h"
#include "subs_matrix.h"


class Alignment 
{

public:
	Alignment(size_t size, const SubstitutionMatrix& sm);

	typedef Matrix<AlnScoreType> ScoreMatrix;

	AlnScoreType globalAlignment(const std::string& consensus,
								 const std::vector<std::string>& reads);

	AlnScoreType addDeletion(unsigned int letterIndex) const;
	AlnScoreType addSubstitution(unsigned int letterIndex,
						   		 char base, const std::vector<std::string>& reads) const;
	AlnScoreType addInsertion(unsigned int positionIndex,
						   	  char base, const std::vector<std::string>& reads) const;

private:
	std::vector<ScoreMatrix> _forwardScores;
	std::vector<ScoreMatrix> _reverseScores;
	const SubstitutionMatrix& _subsMatrix;

	AlnScoreType getScoringMatrix(const std::string& v, const std::string& w,
							      ScoreMatrix& scoreMat);
};
