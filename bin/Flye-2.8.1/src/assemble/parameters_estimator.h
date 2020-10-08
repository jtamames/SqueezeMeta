//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"
#include <limits>

class ParametersEstimator
{
public:
	ParametersEstimator(const SequenceContainer& seqContainer,
						const VertexIndex& vertexIndex, size_t genomeSize):
		_vertexIndex(vertexIndex), 
		_seqContainer(seqContainer),
		_genomeSize(genomeSize),
		_minKmerCount(std::numeric_limits<size_t>::max())
	{}

	void    estimateMinKmerCount();
	size_t  genomeSizeEstimate();
	size_t 	minKmerCount() {return _minKmerCount;}
private:

	const VertexIndex& _vertexIndex;
	const SequenceContainer& _seqContainer;
	const size_t _genomeSize;
	size_t _takenKmers;
	size_t _minKmerCount;
};
