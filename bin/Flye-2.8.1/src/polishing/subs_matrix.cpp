//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <sstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cassert>

#include "subs_matrix.h"
#include "utility.h"

namespace
{
	static std::vector<size_t> table;
	size_t dnaToId(char c)
	{
		return table[(size_t)c];
	}
	struct TableFiller
	{
		TableFiller()
		{
			static bool tableFilled = false;
			if (!tableFilled)
			{
				tableFilled = true;
				table.assign(256, -1);	//256 chars
				table[(size_t)'A'] = 0;
				table[(size_t)'a'] = 0;
				table[(size_t)'C'] = 1;
				table[(size_t)'c'] = 1;
				table[(size_t)'T'] = 2;
				table[(size_t)'t'] = 2;
				table[(size_t)'G'] = 3;
				table[(size_t)'g'] = 3;
				table[(size_t)'-'] = 4;
			}
		}
	};
	TableFiller filler;


	static const size_t MIN_HOPO = 1;
	static const size_t MAX_HOPO = 20;
	static const size_t NUM_HOPO_STATES = 128;
	static const size_t NUM_HOPO_OBS = 65536;
	static const double MIN_HOPO_PROB = 0.001f;
	static const double ZERO_HOPO_PROB = 0.0000000001f;

	static const int32_t SCORE_MULT = 2 << 16;
	AlnScoreType probToScore(double prob)
	{
		//std::cout << prob << " " << (int32_t)std::round(std::log(prob) * (double)SCORE_MULT) << "\n";
		return std::round(std::log(prob) * (double)SCORE_MULT);
	}
}

SubstitutionMatrix::SubstitutionMatrix(const std::string& path)
{	
	_matrix.assign(MAX_CHAR * MAX_CHAR, AlnScoreType(0));
	this->loadMatrix(path);
}

void SubstitutionMatrix::loadMatrix(const std::string& path) 
{
	std::string line;
	std::ifstream file(path);

	if (!file.is_open()) 
	{
		throw std::runtime_error("Can't open substitution matrix");
	}
	while (std::getline(file, line))
	{
		if (line.empty()) continue;
		if (line[line.length() - 1] == '\r') line.pop_back();

		auto items = splitString(line, '\t');

		if (items[0] == "mat") 
		{
			char nuc = items[1][0];
			double probability = std::stod(items[2]);
			this->setScore(nuc, nuc, probToScore(probability));
		}
		else if (items[0] == "mis") 
		{
			char x = items[1][0];
			char y = items[1][items[1].size() - 1];
			double probability = std::stod(items[2]);
			this->setScore(x, y, probToScore(probability));
		}
		else if (items[0] == "del") 
		{
			char x = items[1][0];
			double probability = std::stod(items[2]);
			this->setScore(x, '-', probToScore(probability));
		}
		else if (items[0] == "ins") 
		{
			char y = items[1][0];
			double probability = std::stod(items[2]);
			this->setScore('-', y, probToScore(probability));
		}	
	}
}


namespace
{
	//converts condensed hopo string to a normal nucleotide string
	//e.g. 10A5C3T to AAAAAAAAAACCCCCTTT
	std::string expandHopo(const std::string& hopo)
	{
		std::string buf;
		std::string result;
		for (char c : hopo)
		{
			if (!std::isalpha(c))
			{
				buf += c;
			}
			else
			{
				int runLen = std::stoi(buf);
				result += std::string(runLen, c);
				buf.clear();
			}
		}
		return result;
	}
}

//constructs state for a given nucleotide and length (e.g. 5C or 3A)
HopoMatrix::State::State(char nucl, uint32_t length):
	nucl(nucl), length(length)
{
	 id = std::min(length, (uint32_t)MAX_HOPO) + (MAX_HOPO + 1) * dnaToId(nucl);
}

		
//constructs State from dna homopolymer string
HopoMatrix::State::State(const std::string& str, size_t start, size_t end)
{
	if (end == std::string::npos) end = str.length();
	assert(end - start > 0);

	//std::cerr << str.substr(start, end - start) << std::endl;
	size_t runLength = 0;
	char runNucl = -1;
	for (size_t i = start; i < end; ++i)
	{
		if (str[i] != '-')
		{
			if (runNucl == -1) runNucl = str[i];
			if (str[i] != runNucl) throw std::runtime_error("Wrong homopolymer");
			++runLength;
		}
	}
	if (runLength == 0) throw std::runtime_error("Wrong homopolymer");

	id = std::min(runLength, MAX_HOPO) + (MAX_HOPO + 1) * dnaToId(runNucl);
	nucl = runNucl;
	length = runLength;
}

//converts nucleotide string into observation, given the 
//main homopolymer nucleotide, e.g. AAACAAGA is converted to
//6A2X
HopoMatrix::Observation HopoMatrix::strToObs(char mainNucl, 
											 const std::string& dnaStr, 
											 size_t begin, size_t end)
{
	size_t counts[] = {0, 0};
	uint32_t result = 0;

	if (end == std::string::npos) end = dnaStr.length();
	for (size_t pos = begin; pos < end; ++pos)
	{
		if (dnaStr[pos] != '-')
		{
			++counts[dnaStr[pos] == mainNucl];
		}
	}

	for (size_t i = 0; i < 2; ++i)
	{
		counts[i] = std::min(counts[i], MAX_HOPO);
		result <<= 4;
		result += (uint32_t)counts[i];
	}
	return Observation(result);
}

/*
std::string HopoMatrix::obsToStr(HopoMatrix::Observation obs)
{
	//ACTG
	char NUCLS[] = "GTCA";
	std::string result;
	for (size_t i = 0; i < 4; ++i)
	{
		int num = obs.id & MAX_RUN;
		obs.id >>= 4;
		if (num) result += std::to_string(num) + NUCLS[i];
	}
	return result;
}*/

HopoMatrix::HopoMatrix(const std::string& fileName)
{
	for (size_t i = 0; i < NUM_HOPO_STATES; ++i)
	{
		_observationProbs.emplace_back(NUM_HOPO_OBS, probToScore(MIN_HOPO_PROB));
	}
	_genomeProbs.assign(NUM_HOPO_STATES, probToScore(MIN_HOPO_PROB));
	this->loadMatrix(fileName);
}


//returns all observations from the training set for a given state
HopoMatrix::ObsVector HopoMatrix::knownObservations(State state) const
{
	ObsVector obsVector;
	for (uint32_t obsId = 0; obsId < NUM_HOPO_OBS; ++obsId)
	{
		if (_observationProbs[state.id][obsId] > probToScore(MIN_HOPO_PROB))
		{
			obsVector.emplace_back(obsId);
		}
	}
	return obsVector;
}

//loads homopolymer matrix from .mat file
void HopoMatrix::loadMatrix(const std::string& fileName)
{
	std::ifstream fin(fileName);
	std::string buffer;
	if (!fin.is_open()) 
	{
		throw std::runtime_error("Can't open homopolymer matrix");
	}

	std::vector<size_t> nucleotideFreq(5, 0);
	std::vector<std::vector<int>> observationsFreq;
	for (size_t i = 0; i < NUM_HOPO_STATES; ++i)
	{
		observationsFreq.push_back(std::vector<int>(NUM_HOPO_OBS, 0));
	}

	while (std::getline(fin, buffer))
	{
		if (buffer.empty()) continue;
		if (buffer[buffer.length() - 1] == '\r') buffer.pop_back();

		auto tokens = splitString(buffer, '\t');
		tokens[0].pop_back();
		State state = State(expandHopo(tokens[0]));
		for (size_t i = 1; i < tokens.size(); ++i)
		{
			auto obsTokens = splitString(tokens[i], '=');
			Observation obs = strToObs(state.nucl, expandHopo(obsTokens[0]));
			observationsFreq[state.id][obs.id] += std::stoi(obsTokens[1]);
			nucleotideFreq[dnaToId(state.nucl)] += std::stoi(obsTokens[1]);
		}
	}


	for (char nucl : std::string("ACGT"))
	{
		for (size_t runLen = MIN_HOPO; runLen <= MAX_HOPO; ++runLen)
		{
			State state(nucl, runLen);
			int sumFreq = 0;
			for (size_t j = 0; j < NUM_HOPO_OBS; ++j)
			{
				sumFreq += observationsFreq[state.id][j];
			}
			double prob = (double)sumFreq / nucleotideFreq[dnaToId(state.nucl)];
			_genomeProbs[state.id] = probToScore(std::max(prob, ZERO_HOPO_PROB));
			//std::cerr << state.length << state.nucl << "\t" << _genomeProbs[state.id] 
			//		  << "\t" << sumFreq << "\t" << nucleotideFreq[dnaToId(state.nucl)]
			//		  << std::endl;

			if (sumFreq == 0) continue;
			for (size_t j = 0; j < NUM_HOPO_OBS; ++j)
			{
				double prob = (double)observationsFreq[state.id][j] / sumFreq;
				_observationProbs[state.id][j] = 
									probToScore(std::max(prob, MIN_HOPO_PROB));
			}
		}
	}
}
