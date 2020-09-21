//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "alignment.h"
#include <chrono>


Alignment::Alignment(size_t size, const SubstitutionMatrix& sm):
	_forwardScores(size),
	_reverseScores(size),
	_subsMatrix(sm)
{ 
}


AlnScoreType Alignment::globalAlignment(const std::string& consensus,
							 			const std::vector<std::string>& reads)
{
	AlnScoreType finalScore = 0;
	for (size_t readId = 0; readId < _forwardScores.size(); ++readId)
	{
		unsigned int x = consensus.size() + 1;
		unsigned int y = reads[readId].size() + 1;
		ScoreMatrix scoreMat(x, y, 0);
			
		AlnScoreType score = this->getScoringMatrix(consensus, reads[readId], 
													  scoreMat);
		_forwardScores[readId] = std::move(scoreMat);

		//The reverse alignment is similar, but we need
		//the scoring matrix with of the reverse alignment (I guess)
		std::string revConsensus(consensus.rbegin(), consensus.rend());
		std::string revRead(reads[readId].rbegin(), reads[readId].rend());

		ScoreMatrix scoreMatRev(x, y, 0);
		this->getScoringMatrix(revConsensus, revRead, scoreMatRev);
		_reverseScores[readId] = std::move(scoreMatRev);

		finalScore += score;
	}
	return finalScore;
}

AlnScoreType Alignment::addDeletion(unsigned int letterIndex) const
{
	AlnScoreType finalScore = 0;
	for (size_t readId = 0; readId < _forwardScores.size(); ++readId)
	{
		const ScoreMatrix& forwardScore = _forwardScores[readId];
		const ScoreMatrix& reverseScore = _reverseScores[readId];

		//Note: We subtract 2 because of zero indexing and an extra added row and column count
		//unsigned int index = (reverseScore.nrows() - 1) - letterIndex;
		size_t frontRow = letterIndex - 1;
		size_t revRow = reverseScore.nrows() - 1 - letterIndex;
		
		AlnScoreType maxVal = std::numeric_limits<AlnScoreType>::lowest();
		for (size_t col = 0; col < forwardScore.ncols(); ++col)
		{
			size_t backCol = forwardScore.ncols() - col - 1;
			AlnScoreType sum = forwardScore.at(frontRow, col) + 
							reverseScore.at(revRow, backCol);
			maxVal = std::max(maxVal, sum);
		}
		finalScore += maxVal;
	}
	return finalScore;
}


//AlnScoreType Alignment::addSubstitution(unsigned int wordIndex, 
//							   			unsigned int letterIndex,
//							   			char base, const std::string& read) 

AlnScoreType Alignment::addSubstitution(unsigned int letterIndex, char base, 
										const std::vector<std::string>& reads) const
{
	AlnScoreType finalScore = 0;
	for (size_t readId = 0; readId < reads.size(); ++readId)
	{
		//LetterIndex must start with 1 and go until (row.size - 1)
		const ScoreMatrix& forwardScore = _forwardScores[readId];
		const ScoreMatrix& reverseScore = _reverseScores[readId];

		size_t frontRow = letterIndex - 1;
		size_t revRow = reverseScore.nrows() - 1 - letterIndex;

		std::vector<AlnScoreType> sub(reads[readId].size() + 1);
		sub[0] = forwardScore.at(frontRow, 0) + _subsMatrix.getScore(base, '-');
		for (size_t i = 0; i < reads[readId].size(); ++i)
		{
			AlnScoreType match = forwardScore.at(frontRow, i) + 
							_subsMatrix.getScore(base, reads[readId][i]);
			AlnScoreType ins = forwardScore.at(frontRow, i + 1) + 
							_subsMatrix.getScore(base, '-');
			sub[i + 1] = std::max(match, ins);
		}

		AlnScoreType maxVal = std::numeric_limits<AlnScoreType>::lowest();
		for (size_t col = 0; col < forwardScore.ncols(); ++col)
		{
			size_t backCol = forwardScore.ncols() - col - 1;
			AlnScoreType sum = sub[col] + reverseScore.at(revRow, backCol);
			maxVal = std::max(maxVal, sum);
		}
		finalScore += maxVal;
	}
	return finalScore;
}


AlnScoreType Alignment::addInsertion(unsigned int pos, char base, 
									 const std::vector<std::string>& reads) const
{
	AlnScoreType finalScore = 0;
	for (size_t readId = 0; readId < reads.size(); ++readId)
	{
		//LetterIndex must start with 1 and go until (row.size - 1)
		const ScoreMatrix& forwardScore = _forwardScores[readId];
		const ScoreMatrix& reverseScore = _reverseScores[readId];

		size_t frontRow = pos - 1;
		size_t revRow = reverseScore.nrows() - pos;

		std::vector<AlnScoreType> sub(reads[readId].size() + 1);
		sub[0] = forwardScore.at(frontRow, 0) + _subsMatrix.getScore(base, '-');
		for (size_t i = 0; i < reads[readId].size(); ++i)
		{
			AlnScoreType match = forwardScore.at(frontRow, i) + 
							_subsMatrix.getScore(base, reads[readId][i]);
			AlnScoreType ins = forwardScore.at(frontRow, i + 1) + 
							_subsMatrix.getScore(base, '-');
			sub[i + 1] = std::max(match, ins);
		}

		AlnScoreType maxVal = std::numeric_limits<AlnScoreType>::lowest();
		for (size_t col = 0; col < forwardScore.ncols(); ++col)
		{
			size_t backCol = forwardScore.ncols() - col - 1;
			AlnScoreType sum = sub[col] + reverseScore.at(revRow, backCol);
			maxVal = std::max(maxVal, sum);
		}
		finalScore += maxVal;
	}
	return finalScore;
}


AlnScoreType Alignment::getScoringMatrix(const std::string& v, 
										 const std::string& w,
								  		 ScoreMatrix& scoreMat) 
{
	AlnScoreType score = 0;
	
	for (size_t i = 0; i < v.size(); i++) 
	{
		AlnScoreType score = _subsMatrix.getScore(v[i], '-');
		scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + score;
	}


	for (size_t i = 0; i < w.size(); i++) {
		AlnScoreType score = _subsMatrix.getScore('-', w[i]);
		scoreMat.at(0, i + 1) = scoreMat.at(0, i) + score;
	}


	for (size_t i = 1; i < v.size() + 1; i++)
	{
		char key1 = v[i - 1];
		for (size_t j = 1; j < w.size() + 1; j++) 
		{
			char key2 = w[j - 1];

			AlnScoreType left = scoreMat.at(i, j - 1) + 
							_subsMatrix.getScore('-', key2);
			AlnScoreType up = scoreMat.at(i - 1, j) + 
							_subsMatrix.getScore(key1, '-');
			score = std::max(left, up);

			AlnScoreType cross = scoreMat.at(i - 1, j - 1) + 
							_subsMatrix.getScore(key1, key2);
			score = std::max(score, cross);
			scoreMat.at(i, j) = score;
		}
	}

	return score;
}
