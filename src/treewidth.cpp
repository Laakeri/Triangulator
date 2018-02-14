#include "treewidth.hpp"
#include "graph.hpp"
#include "MCS.hpp"
#include "global.hpp"
#include "PMC.hpp"

#include <vector>
#include <cassert>
#include <algorithm>
#include <queue>
#include <tuple>
#include <iostream>

#define F first
#define S second

using namespace std;

namespace TreeWidth {

namespace SolveAtom {

pair<int, vector<pair<int, int> > > solveAtom1(const Graph& G, int lowerBound, int okWidth, int upperBound) {
	assert(lowerBound >= 0);
	assert(okWidth <= upperBound);
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	int fillWidth = MCS::getTreeWidth(gFill);
	assert(fillWidth >= lowerBound);
	if (fillWidth <= okWidth) return {fillWidth, fill.F};
	if (fill.F.size() <= 1) {
		if (fillWidth > upperBound) return {-1, FAILSOLUTION};
		return {fillWidth, fill.F};
	}
	auto fillO = PMC::solve(G, (int64_t)lowerBound, (int64_t)upperBound, PROBLEM_TREEWIDTH);
	if (fillO == FAILSOLUTION) return {-1, FAILSOLUTION};
	Graph gFillO = G;
	gFillO.addEdges(fillO);
	int width = MCS::getTreeWidth(gFillO);
	assert(width >= lowerBound);
	assert(width <= upperBound);
	return {width, fillO};
}
vector<pair<int, int> > greedyDegree2(Graph& G) {
	queue<int> d2;
	for (int i = 0; i < G.n; i++) {
		assert(G.g[i].size() > 1);
		if (G.g[i].size() == 2) {
			d2.push(i);
		}
	}
	vector<pair<int, int> > sol;
	while (!d2.empty()) {
		int x = d2.front();
		d2.pop();
		if (G.g[x].size() != 2) continue;
		int n1 = G.g[x][0];
		int n2 = G.g[x][1];
		if (!G.hasEdge(n1, n2)) {
			sol.push_back({n1, n2});
			G.addEdge(n1, n2);
		}
	}
	return sol;
}
pair<int, int> nbClique(Graph& G) {
	for (int x = 0; x < G.n; x++) {
		assert(G.g[x].size() > 2);
		pair<int, int> fo = {-1, -1};
		bool fail = false;
		for (int i = 0; i < (int)G.g[x].size(); i++) {
			for (int ii = i+1; ii < (int)G.g[x].size(); ii++) {
				if (!G.hasEdge(G.g[x][i], G.g[x][ii])) {
					if (fo.F != -1) {
						fail = true;
						break;
					}
					else {
						fo = {G.g[x][i], G.g[x][ii]};
					}
				}
			}
			if (fail) break;
		}
		if (!fail) {
			assert(fo.F != -1);
			return fo;
		}
	}
	return {-1, -1};
}
pair<int, vector<pair<int, int> > > solveAtom0(Graph& G, int lowerBound, int okWidth, int upperBound) {
	assert(lowerBound >= 0);
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	if ((int64_t)G.n*(int64_t)(G.n-1) == (int64_t)G.m*2) {
		if (G.n-1 <= upperBound) {
			return {G.n-1, {}};
		}
		else {
			return {-1, FAILSOLUTION};
		}
	}
	if (G.hasCycle()) {
		lowerBound = max(lowerBound, 2);
	}
	else {
		upperBound = min(upperBound, 1);
	}
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);

	if (G.n == 3) {
		if (G.m == 3) {
			return {2, {}};
		}
		else {
			return {1, {}};
		}
	}
	bool hasd2 = false;
	for (int i = 0; i < G.n; i++) {
		assert(G.g[i].size() > 1);
		if (G.g[i].size() == 2) {
			hasd2 = true;
			break;
		}
	}
	if (hasd2) {
		auto sol0 = greedyDegree2(G);
		auto sol = SolveGraph::solveGraph(G, lowerBound, okWidth, upperBound);
		if (sol.S == FAILSOLUTION) return sol;
		sol.S.insert(sol.S.end(), sol0.begin(), sol0.end());
		assert(sol.F <= upperBound);
		return sol;
	}
	auto nbC = nbClique(G);
	if (nbC.F != -1) {
		auto sol0 = nbC;
		G.addEdge(sol0.F, sol0.S);
		auto sol = SolveGraph::solveGraph(G, lowerBound, okWidth, upperBound);
		if (sol.S == FAILSOLUTION) return sol;
		sol.S.push_back(sol0);
		assert(sol.F <= upperBound);
		return sol;
	}
	auto ret = solveAtom1(G, lowerBound, okWidth, upperBound);
	if (ret.S != FAILSOLUTION) assert(ret.F <= upperBound);
	return ret;
}

}

namespace SolveGraph {

// Find a fill-in with width at most okWidth or the best solution that is at most upperBound. There is no solution better than lowerBound.
pair<int, vector<pair<int, int> > > solveGraph(const Graph& G, int lowerBound, int okWidth, int upperBound) {
	assert(lowerBound >= 0);
	assert(okWidth <= upperBound);
	if (G.hasCycle()) {
		lowerBound = max(lowerBound, 2);
	}
	else {
		upperBound = min(upperBound, 1);
	}
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	int fillWidth = MCS::getTreeWidth(gFill);
	assert(fillWidth >= lowerBound);
	if (fillWidth <= okWidth) return {fillWidth, fill.F};
	
	if (fill.F.size() <= 1) {
		if (fillWidth > upperBound) return {-1, FAILSOLUTION};
		return {fillWidth, fill.F};
	}
	
	vector<Graph> atoms;
	int maxCliqueSeparator;
	tie(maxCliqueSeparator, atoms) = MCS::getAtoms(G, gFill, fill);
	lowerBound = max(lowerBound, maxCliqueSeparator-1);
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	if (fillWidth <= okWidth) return {fillWidth, fill.F};
	
	auto cmp = [](const Graph& a, const Graph& b) {
		return a.n > b.n;
	};
	sort(atoms.begin(), atoms.end(), cmp);
	pair<int, vector<pair<int, int> > > sol = {lowerBound, {}};
	for (int i = 0; i < (int)atoms.size(); i++) {
		pair<int, vector<pair<int, int> > > tSol;
		if (i+1 < (int)atoms.size()) {
			tSol = SolveAtom::solveAtom0(atoms[i], 0, sol.F, min(upperBound, fillWidth-1));
		}
		else {
			tSol = SolveAtom::solveAtom0(atoms[i], 0, sol.F, min(upperBound, fillWidth-1));
		}
		if (tSol.S == FAILSOLUTION) {
			if (fillWidth <= upperBound) {
				return {fillWidth, fill.F};
			}
			else {
				return {-1, FAILSOLUTION};
			}
		}
		sol.F = max(sol.F, tSol.F);
		tSol.S = atoms[i].mapBack(tSol.S);
		sol.S.insert(sol.S.end(), tSol.S.begin(), tSol.S.end());
	}
	assert(sol.F <= upperBound);
	return sol;
}

}


}