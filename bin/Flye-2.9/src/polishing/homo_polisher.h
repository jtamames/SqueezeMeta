//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "subs_matrix.h"
#include "bubble.h"

class HomoPolisher
{
public:
	HomoPolisher(const SubstitutionMatrix& subsMatrix,
				 const HopoMatrix& hopoMatrix):
		_subsMatrix(subsMatrix), _hopoMatrix(hopoMatrix)
	{}
	void polishBubble(Bubble& bubble) const;

private:
	size_t mostLikelyLen(char nucleotide, 
						 const HopoMatrix::ObsVector& obs) const;
	size_t compareTopTwo(char nucleotide, size_t firstChoice, 
						 size_t secondChoice,
				   		 const HopoMatrix::ObsVector& observations) const;
	AlnScoreType likelihood(HopoMatrix::State state, 
					  		const HopoMatrix::ObsVector& observations) const;

	const SubstitutionMatrix& _subsMatrix;
	const HopoMatrix& _hopoMatrix;
};
