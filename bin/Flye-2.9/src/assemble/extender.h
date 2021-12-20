//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <deque>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../sequence/consensus_generator.h"
#include "chimera.h"

class Extender
{
public:
	Extender(const SequenceContainer& readsContainer, 
			 OverlapContainer& ovlpContainer,
			 int safeOverlap):
		_safeOverlap(safeOverlap),
		_readsContainer(readsContainer), 
		_ovlpContainer(ovlpContainer),
		_chimDetector(readsContainer, ovlpContainer)
	{}

	void assembleDisjointigs();
	const std::vector<ContigPath>& getDisjointigPaths() const
		{return _disjointigPaths;}

private:
	struct ExtensionInfo
	{
		ExtensionInfo(): leftTip(false), rightTip(false),
			numSuspicious(0), meanOverlaps(0), stepsToTurn(0),
			assembledLength(0), singleton(false),
			avgOverlapSize(0), minOverlapSize(0),
			//leftAsmOverlap(0), rightAsmOverlap(0),
			shortExtensions(0) {}

		std::vector<FastaRecord::Id> reads;
		bool leftTip;
		bool rightTip;
		int  numSuspicious;
		int  meanOverlaps;
		int  stepsToTurn;
		int  assembledLength;
		bool singleton;
		int  avgOverlapSize;
		int  minOverlapSize;
		//int  leftAsmOverlap;
		//int  rightAsmOverlap;
		int  shortExtensions;
	};

	const int _safeOverlap;

	ExtensionInfo extendDisjointig(FastaRecord::Id startingRead);
	//int   countRightExtensions(FastaRecord::Id readId) const;
	int   countRightExtensions(const std::vector<OverlapRange>&) const;
	int   countLeftExtensions(const std::vector<OverlapRange>&) const;
	bool  extendsRight(const OverlapRange& ovlp) const;
	bool  extendsLeft(const OverlapRange& ovlp) const;
	void  convertToDisjointigs();
	std::vector<FastaRecord::Id> 
		getInnerReads(const std::vector<OverlapRange>& ovlps);

	const SequenceContainer& _readsContainer;
	OverlapContainer& _ovlpContainer;
	ChimeraDetector   _chimDetector;

	std::vector<ExtensionInfo> 	_readLists;
	std::vector<ContigPath> 	_disjointigPaths;
	cuckoohash_map<FastaRecord::Id, size_t>  	_innerReads;
};
