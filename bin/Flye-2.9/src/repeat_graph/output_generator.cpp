//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "output_generator.h"
#include "../sequence/consensus_generator.h"
#include <iomanip>


//Generates FASTA from the given graph paths
std::vector<FastaRecord> OutputGenerator::
	generatePathSequences(const std::vector<UnbranchingPath>& paths) const
{
	std::vector<FastaRecord> contigSequences;

	for (auto& contig : paths)
	{
		//As each edge might correspond to multiple sequences,
		//we need to select them so as to minimize the
		//number of original contigs (that were used to build the graph)
		std::unordered_map<FastaRecord::Id, int> seqIdFreq;
		for (auto& edge : contig.path) 
		{
			std::unordered_set<FastaRecord::Id> edgeSeqIds;
			for (auto& seg: edge->seqSegments) 
			{
				edgeSeqIds.insert(seg.origSeqId);
			}
			for (auto& seqId : edgeSeqIds)
			{
				seqIdFreq[seqId] += 1;
			}
		}

		std::string nucSequence;
		//ContigPath contigPath;
		//contigPath.name = contig.name();
		//int32_t prevFlank = 0;
		//int32_t prevSubLength = 0;

		for (size_t i = 0; i < contig.path.size(); ++i) 
		{
			if (contig.path[i]->seqSegments.empty()) 
			{
				throw std::runtime_error("Edge without sequence");
			}

			//get the sequence with maximum frequency
			EdgeSequence* bestSegment = nullptr;
			for (auto& seg : contig.path[i]->seqSegments)
			{
				if (!bestSegment || 
					seqIdFreq[seg.origSeqId] > seqIdFreq[bestSegment->origSeqId])
				{
					bestSegment = &seg;
				}
			}
			if (bestSegment->seqLen == 0) continue;
			nucSequence += _graph.edgeSequences()
								.getSeq(bestSegment->edgeSeqId).str();
		}
		//contigParts.push_back(contigPath);
		contigSequences.push_back({DnaSequence(nucSequence), contig.name(), 
								  FastaRecord::ID_NONE});
	}

	return contigSequences;
}

void OutputGenerator::outputFasta(const std::vector<UnbranchingPath>& paths,
								  const std::string& filename)
{
	std::vector<UnbranchingPath> posStrandPaths;
	for (auto& path : paths)
	{
		if (path.id.strand()) posStrandPaths.push_back(path);
	}
	SequenceContainer::writeFasta(this->generatePathSequences(posStrandPaths), 
								  filename);
}

void OutputGenerator::outputGfa(const std::vector<UnbranchingPath>& paths,
							    const std::string& filename)
{
	auto sequences = this->generatePathSequences(paths);
	std::unordered_map<GraphEdge*, const UnbranchingPath*> edgeToPath;
	for (auto& path : paths)
	{
		edgeToPath[path.path.front()] = &path;
	}

	auto edgeConnections = _aligner.getEdgeConnectivity();

	Logger::get().debug() << "Writing Gfa";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);

	fprintf(fout, "H\tVN:Z:1.0\n");
	for (size_t i = 0; i < paths.size(); ++i)
	{
		if (!paths[i].id.strand()) continue;

		//size_t kmerCount = sequences[i].sequence.length() * paths[i].meanCoverage;
		fprintf(fout, "S\t%s\t%s\tdp:i:%d\n", paths[i].name().c_str(), 
				sequences[i].sequence.str().c_str(), (int)paths[i].meanCoverage);
	}

	//make sure that if there are nodes with one incoming and one outgoing
	//edge, they are connected. Most relevant to the circular contigs
	for (auto& node : _graph.iterNodes())
	{
		if (node->inEdges.size() == 1 && node->outEdges.size() == 1)
		{
			//initialize to zero
			edgeConnections[node->inEdges.front()][node->outEdges.front()];
		}
	}

	std::unordered_set<std::pair<GraphEdge*, GraphEdge*>, pairhash> usedPairs;
	for (size_t i = 0; i < paths.size(); ++i)
	{
		for (auto& outEdgeIt : edgeConnections[paths[i].path.back()])
		{
			auto outPath = edgeToPath[outEdgeIt.first];
			GraphEdge* edgeLeft = paths[i].path.back();
			GraphEdge* edgeRight = outEdgeIt.first;

			if (usedPairs.count(std::make_pair(edgeLeft, edgeRight))) continue;
			usedPairs.insert(std::make_pair(edgeLeft, edgeRight));
			usedPairs.insert(std::make_pair(_graph.complementEdge(edgeRight), 
											_graph.complementEdge(edgeLeft)));

			std::string leftSign = paths[i].id.strand() ? "+" :"-";
			std::string leftName = paths[i].nameUnsigned();

			std::string rightSign = outPath->id.strand() ? "+" :"-";
			std::string rightName = outPath->nameUnsigned();

			fprintf(fout, "L\t%s\t%s\t%s\t%s\t0M\tRC:i:%d\n", leftName.c_str(), 
					leftSign.c_str(), rightName.c_str(), rightSign.c_str(), outEdgeIt.second);
		}
	}
	fclose(fout);
}

/*void OutputGenerator::outputGfaCompact(const std::vector<UnbranchingPath>& paths,
							    	   const std::string& filename)
{
	auto sequences = this->generatePathSequences(paths);
	std::unordered_map<GraphNode*, std::vector<const UnbranchingPath*>> leftNodes;
	std::unordered_map<GraphNode*, std::vector<const UnbranchingPath*>> rightNodes;
	for (auto& path : paths)
	{
		leftNodes[path.path.front()->nodeLeft].push_back(&path);
		rightNodes[path.path.back()->nodeRight].push_back(&path);
	}

	Logger::get().debug() << "Writing Gfa";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);

	fprintf(fout, "H\tVN:Z:1.0\n");
	for (size_t i = 0; i < paths.size(); ++i)
	{
		if (!paths[i].id.strand()) continue;

		//size_t kmerCount = sequences[i].sequence.length() * paths[i].meanCoverage;
		fprintf(fout, "S\t%s\t%s\tdp:i:%d\n", paths[i].name().c_str(), 
				sequences[i].sequence.str().c_str(), (int)paths[i].meanCoverage);
	}

	std::unordered_map<GraphNode*, int> signedNodeIds;
	for (auto& nodeIt : leftNodes)
	{
		if (rightNodes.count(nodeIt.first))
		{
			if (!signedNodeIds.count(nodeIt.first))
			{
				signedNodeIds[nodeIt.first] = (int)nodeIt.first->nodeId;
				GraphNode* complNode = _graph.complementNode(nodeIt.first);
				if (complNode != nodeIt.first) 
				{
					signedNodeIds[complNode] = -(int)nodeIt.first->nodeId;
				}

				std::string nodeId = "node_" + std::to_string(nodeIt.first->nodeId);
				fprintf(fout, "S\t%s\t*\n", nodeId.c_str());
			}
		}
	}

	for (auto& path : paths)
	{
		std::string ctgSign = path.id.strand() ? "+" :"-";
		std::string ctgName = path.nameUnsigned();

		GraphNode* nodeLeft = path.path.front()->nodeLeft;
		GraphNode* nodeRight = path.path.back()->nodeRight;

		if (rightNodes.count(nodeLeft))
		{
			std::string nodeId = "node_" + std::to_string(abs(signedNodeIds.at(nodeLeft)));
			std::string nodeSign = signedNodeIds.at(nodeLeft) > 0 ? "+" : "-";
			fprintf(fout, "L\t%s\t%s\t%s\t%s\t0M\n", nodeId.c_str(), 
					nodeSign.c_str(), ctgName.c_str(), ctgSign.c_str());
		}

		if (leftNodes.count(nodeRight))
		{
			std::string nodeId = "node_" + std::to_string(abs(signedNodeIds.at(nodeRight)));
			std::string nodeSign = signedNodeIds.at(nodeRight) > 0 ? "+" : "-";
			fprintf(fout, "L\t%s\t%s\t%s\t%s\t0M\n", ctgName.c_str(), 
					ctgSign.c_str(), nodeId.c_str(), nodeSign.c_str());
		}
	}
}*/

void OutputGenerator::outputDot(const std::vector<UnbranchingPath>& paths,
								const std::string& filename)
{
	Logger::get().debug() << "Writing Dot";

	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);

	fout << "digraph {\n";
	fout << "nodesep = 0.5;\n";
	fout << "node [shape = circle, label = \"\", height = 0.3];\n";
	
	///re-enumerating helper functions
	std::unordered_map<GraphNode*, int> nodeIds;
	int nextNodeId = 0;
	auto nodeToId = [&nodeIds, &nextNodeId](GraphNode* node)
	{
		if (!nodeIds.count(node))
		{
			nodeIds[node] = nextNodeId++;
		}
		return nodeIds[node];
	};

	for (auto& node : _graph.iterNodes())
	{
		if (node->isTelomere())
		{
			fout << "\"" << nodeToId(node) 
				<< "\" [style = \"filled\", fillcolor = \"grey\"];\n";
		}
	}

	//coloring repeat clusters
	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue1", "darkorchid", "aquamarine1", 
								  "darkgoldenrod1", "deepskyblue1", 
								  "darkolivegreen3"};
	std::vector<GraphNode*> dfsStack;
	std::unordered_set<GraphNode*> visited;
	std::unordered_map<GraphEdge*, std::string> edgeColors;
	size_t colorId = 0;
	for (auto& node : _graph.iterNodes())
	{
		if (visited.count(node)) continue;
		dfsStack.push_back(node);
		while(!dfsStack.empty())
		{
			auto curNode = dfsStack.back(); 
			dfsStack.pop_back();
			if (visited.count(curNode)) continue;

			visited.insert(curNode);
			visited.insert(_graph.complementNode(curNode));
			
			for (auto adjEdge: curNode->outEdges)
			{
				if (adjEdge->isRepetitive())
				{
					edgeColors[adjEdge] = COLORS[colorId];
					edgeColors[_graph.complementEdge(adjEdge)] = COLORS[colorId];
					if (!visited.count(adjEdge->nodeRight))
					{
						dfsStack.push_back(adjEdge->nodeRight);
					}
				}
			}
			for (auto adjEdge: curNode->inEdges)
			{
				if (adjEdge->isRepetitive())
				{
					edgeColors[adjEdge] = COLORS[colorId];
					edgeColors[_graph.complementEdge(adjEdge)] = COLORS[colorId];
					if (!visited.count(adjEdge->nodeLeft))
					{
						dfsStack.push_back(adjEdge->nodeLeft);
					}
				}
			}
		}
		colorId = (colorId + 1) % (sizeof(COLORS) / sizeof(COLORS[0]));
	}

	for (auto& contig : paths)
	{
		std::stringstream lengthStr;
		if (contig.length < 5000)
		{
			lengthStr << std::fixed << std::setprecision(1) 
				<< (float)contig.length / 1000 << "k";
		}
		else
		{
			lengthStr << contig.length / 1000 << "k";
		}
		lengthStr << " " << contig.meanCoverage << "x";

		if (contig.repetitive)
		{
			std::string color = edgeColors[contig.path.front()];
			std::string direction = contig.path.front()->selfComplement ?
									", dir = both" : "";
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft) 
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId() << 
				 "\\l" << lengthStr.str() << "\", color = \"" 
				 << color << "\" " << ", penwidth = 3" << direction << "] ;\n";
		}
		else
		{
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft)
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId()
				 << "\\l" << lengthStr.str() << "\", color = \"black\"] ;\n";
		}
	}

	fout << "}\n";
}
