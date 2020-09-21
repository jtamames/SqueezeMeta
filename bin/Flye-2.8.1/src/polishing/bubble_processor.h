//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <mutex>
#include <fstream>

#include "subs_matrix.h"
#include "bubble.h"
#include "general_polisher.h"
#include "homo_polisher.h"
#include "utility.h"
#include "../common/progress_bar.h"
#include "dinucleotide_fixer.h"


class BubbleProcessor 
{
public:
	BubbleProcessor(const std::string& subsMatPath,
					const std::string& hopoMatrixPath,
					bool  showProgress);
	void polishAll(const std::string& inBubbles, const std::string& outConsensus,
				   int numThreads);
	void enableVerboseOutput(const std::string& filename);

private:
	void parallelWorker();
	void cacheBubbles(int numBubbles);
	void writeBubbles(const std::vector<Bubble>& bubbles);
	void writeLog(const std::vector<Bubble>& bubbles);

	const int BUBBLES_CACHE = 100;

	const SubstitutionMatrix  _subsMatrix;
	const HopoMatrix 		  _hopoMatrix;
	const GeneralPolisher 	  _generalPolisher;
	const HomoPolisher 		  _homoPolisher;
	const DinucleotideFixer	  _dinucFixer;

	ProgressPercent 		  _progress;
	std::mutex				  _stateMutex;
	std::vector<Bubble>		  _cachedBubbles;

	std::ifstream			  _bubblesFile;
	std::ofstream			  _consensusFile;
	std::ofstream			  _logFile;
	bool					  _verbose;
	bool 					  _showProgress;
};
