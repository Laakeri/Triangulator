#ifndef MINFILL_TREEWIDTH_HPP
#define MINFILL_TREEWIDTH_HPP

#include "graph.hpp"

#include <vector>
#include <ostream>

namespace TreeWidth {
	namespace SolveGraph {
		std::pair<int, std::vector<std::pair<int, int> > > solveGraph(const Graph& G, int lowerBound, int okWidth, int upperBound);
	}
}

#endif