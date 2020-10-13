//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <chrono>
#include <thread>
#include <sys/stat.h>

#include "bubble_processor.h"

namespace
{
	size_t fileSize(const std::string& filename)
	{
		struct stat st;
		if (stat(filename.c_str(), &st) != 0) return 0;
		return st.st_size;
	}
}

BubbleProcessor::BubbleProcessor(const std::string& subsMatPath,
								 const std::string& hopoMatrixPath,
								 bool showProgress):
	_subsMatrix(subsMatPath),
	_hopoMatrix(hopoMatrixPath),
	_generalPolisher(_subsMatrix),
	_homoPolisher(_subsMatrix, _hopoMatrix),
	_dinucFixer(_subsMatrix),
	_verbose(false),
	_showProgress(showProgress)
{
}


void BubbleProcessor::polishAll(const std::string& inBubbles, 
								const std::string& outConsensus,
			   					int numThreads)
{
	_cachedBubbles.clear();
	_cachedBubbles.reserve(BUBBLES_CACHE);

	size_t fileLength = fileSize(inBubbles);
	if (!fileLength)
	{
		throw std::runtime_error("Empty bubbles file!");
	}
	_bubblesFile.open(inBubbles);
	if (!_bubblesFile.is_open())
	{
		throw std::runtime_error("Error opening bubbles file");
	}

	_progress.setFinalCount(fileLength);

	_consensusFile.open(outConsensus);
	if (!_consensusFile.is_open())
	{
		throw std::runtime_error("Error opening consensus file");
	}

	std::vector<std::thread> threads(numThreads);
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i] = std::thread(&BubbleProcessor::parallelWorker, this);
	}
	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i].join();
	}
	if (_showProgress) _progress.setDone();
}


void BubbleProcessor::parallelWorker()
{
	const int MAX_BUBBLE = 5000;

	_stateMutex.lock();
	while (true)
	{
		if (_cachedBubbles.empty())
		{
			this->cacheBubbles(BUBBLES_CACHE);
			if(_cachedBubbles.empty())
			{
				_stateMutex.unlock();
				return;
			}
		}

		Bubble bubble = _cachedBubbles.back();
		_cachedBubbles.pop_back();

		if (bubble.candidate.size() < MAX_BUBBLE &&
			bubble.branches.size() > 1)
		{
			_stateMutex.unlock();
			_generalPolisher.polishBubble(bubble);
			//_homoPolisher.polishBubble(bubble);
			_dinucFixer.fixBubble(bubble);
			_stateMutex.lock();
		}
		
		this->writeBubbles({bubble});
		if (_verbose) this->writeLog({bubble});
	}
}


void BubbleProcessor::writeBubbles(const std::vector<Bubble>& bubbles)
{
	for (auto& bubble : bubbles)
	{
		_consensusFile << ">" << bubble.header << " " << bubble.position
			 		   << " " << bubble.branches.size() << std::endl
			 		   << bubble.candidate << std::endl;
	}
}

void BubbleProcessor::enableVerboseOutput(const std::string& filename)
{
	_verbose = true;
	_logFile.open(filename);
	if (!_logFile.is_open())
	{
		throw std::runtime_error("Error opening log file");
	}
}

void BubbleProcessor::writeLog(const std::vector<Bubble>& bubbles)
{
	std::vector<std::string> methods = {"None", "Insertion", "Substitution",
										"Deletion", "Homopolymer"};

	for (auto& bubble : bubbles)
	{
		for (auto& stepInfo : bubble.polishSteps)
		{
			 _logFile << std::fixed
				 << std::setw(22) << std::left << "Consensus: " 
				 << std::right << stepInfo.sequence << std::endl
				 << std::setw(22) << std::left << "Score: " << std::right 
				 << std::setprecision(2) << stepInfo.score << std::endl;

			_logFile << std::endl;
		}
		_logFile << "-----------------\n";
	}
}


void BubbleProcessor::cacheBubbles(int maxRead)
{
	std::string buffer;
	std::string candidate;

	int readBubbles = 0;
	while (!_bubblesFile.eof() && readBubbles < maxRead)
	{
		std::getline(_bubblesFile, buffer);
		if (buffer.empty()) break;

		std::vector<std::string> elems = splitString(buffer, ' ');
		if (elems.size() < 3 || elems[0][0] != '>')
		{
			throw std::runtime_error("Error parsing bubbles file");
		}
		std::getline(_bubblesFile, candidate);
		std::transform(candidate.begin(), candidate.end(), 
				       candidate.begin(), ::toupper);
		
		Bubble bubble;
		bubble.candidate = candidate;
		bubble.header = elems[0].substr(1, std::string::npos);
		bubble.position = std::stoi(elems[1]);
		int numOfReads = std::stoi(elems[2]);

		int count = 0;
		while (count < numOfReads) 
		{
			if (buffer.empty()) break;

			std::getline(_bubblesFile, buffer);
			std::getline(_bubblesFile, buffer);
			std::transform(buffer.begin(), buffer.end(), 
				       	   buffer.begin(), ::toupper);
			bubble.branches.push_back(buffer);
			count++;
		}
		if (count != numOfReads)
		{
			throw std::runtime_error("Error parsing bubbles file");
		}

		_cachedBubbles.push_back(std::move(bubble));
		++readBubbles;
	}

	int64_t filePos = _bubblesFile.tellg();
	if (_showProgress && filePos > 0)
	{
		_progress.setValue(filePos);
	}
}
