//(c) 2016-2019 by Authors
//This file is a part of the Flye program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"




class HaplotypeResolver
{
public:
	HaplotypeResolver(RepeatGraph& graph, ReadAligner& aligner,
					  const SequenceContainer& asmSeqs, 
					  const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), _asmSeqs(asmSeqs), 
		_readSeqs(readSeqs), _nextAltGroupId(2) {}

	void resetEdges();
	int  findHeterozygousLoops();
	int  findHeterozygousBulges();
	int  findRoundabouts();
	int  findSuperbubbles();
	void collapseHaplotypes();

private:
	DnaSequence pathSequence(GraphPath& path);
	void separeteAdjacentEdges(GraphEdge* inEdge, GraphEdge* outEdge);
	void separateDistantEdges(GraphEdge* inEdge, GraphEdge* outEdge,
							  EdgeSequence insSequence, FastaRecord::Id newId);

	struct PathWithScore
	{
		GraphAlignment path;
		int score;
	};

	struct VariantPaths
	{
		VariantPaths(): startEdge(nullptr), endEdge(nullptr) {}
		GraphEdge* startEdge;
		GraphEdge* endEdge;
		std::vector<PathWithScore> altPaths;
		DnaSequence bridgingSequence;
	};

	VariantPaths findVariantSegment(GraphEdge* startEdge, 
									const std::vector<GraphAlignment>& alnignments,
									const std::unordered_set<GraphEdge*>& loopedEdges);

	RepeatGraph& _graph;
	ReadAligner& _aligner;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;

	std::unordered_map<std::pair<GraphEdge*, GraphEdge*>, 
					   DnaSequence, pairhash> _bridgingSeqs;
	int _nextAltGroupId;
};
