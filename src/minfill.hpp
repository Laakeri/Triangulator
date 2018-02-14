#ifndef MINFILL_MINFILL_HPP
#define MINFILL_MINFILL_HPP

#include "graph.hpp"

#include <vector>
#include <ostream>

namespace MinFill {
	namespace SolveGraph {
		std::vector<std::pair<int, int> > solveGraph(const Graph& G, int lowerBound, int upperBound, int depth);
	}
}

#endif