#include "hyper.hpp"
#include "global.hpp"
#include "graph.hpp"
#include "MCS.hpp"
#include "PMC.hpp"
#include "setcover.hpp"

#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <cassert>
#include <set>
#include <algorithm>
#include <iomanip>

#define F first
#define S second

using namespace std;

typedef long double ld;
const ld eps = 1e-9;

namespace Hyper {

namespace SolveAtom {

pair<int, vector<pair<int, int> > > solveAtom(const Graph& G, int lowerBound, int okWidth, int upperBound, const vector<int>& vertMap) {
	assert(lowerBound >= 0);
	assert(okWidth <= upperBound);
	
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	int fillWidth = MCS::getHyperTreeWidth(gFill, vertMap);
	assert(fillWidth >= lowerBound);
	
	if (fillWidth <= okWidth) return {fillWidth, fill.F};
	
	if (fill.F.size() == 0) {
		if (fillWidth <= upperBound) return {fillWidth, fill.F};
		else return {-1, FAILSOLUTION};
	}
	
	PMC::initVertMap(vertMap);
	auto fillO = PMC::solve(G, (int64_t)lowerBound, (int64_t)upperBound, PROBLEM_HYPER);
	if (fillO == FAILSOLUTION) return {-1, FAILSOLUTION};
	Graph gFillO = G;
	gFillO.addEdges(fillO);
	int width = MCS::getHyperTreeWidth(gFillO, vertMap);
	assert(width >= lowerBound);
	assert(width <= upperBound);
	return {width, fillO};
}

pair<ld, vector<pair<int, int> > > solveFracAtom(const Graph& G, ld lowerBound, ld okWidth, ld upperBound, const vector<int>& vertMap) {
	assert(lowerBound > -eps);
	assert(okWidth < upperBound + eps);
	
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	ld fillWidth = MCS::getFracHyperTreeWidth(gFill, vertMap);
	assert(fillWidth > lowerBound - eps);
	
	if (fillWidth < okWidth + eps) return {fillWidth, fill.F};
	
	if (fill.F.size() == 0) {
		if (fillWidth < upperBound + eps) return {fillWidth, fill.F};
		else return {-1, FAILSOLUTION};
	}
	
	PMC::initVertMap(vertMap);
	auto fillO = PMC::solve(G, lowerBound, upperBound, PROBLEM_FRACHYPER);
	if (fillO == FAILSOLUTION) return {-1, FAILSOLUTION};
	Graph gFillO = G;
	gFillO.addEdges(fillO);
	ld width = MCS::getFracHyperTreeWidth(gFillO, vertMap);
	assert(width > lowerBound - eps);
	assert(width < upperBound + eps);
	return {width, fillO};
}

}

namespace SolveGraph {

pair<int, vector<pair<int, int> > > solveGraph(const Graph& G, int lowerBound, int okWidth, int upperBound, const vector<int>& vertMap) {
	assert(lowerBound >= 0);
	assert(okWidth <= upperBound);
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	int fillWidth = MCS::getHyperTreeWidth(gFill, vertMap);
	assert(fillWidth >= lowerBound);
	
	if (fillWidth <= okWidth) return {fillWidth, fill.F};
	
	vector<Graph> atoms = MCS::getAtoms(G, gFill, fill).S;
	
	auto cmp = [](const Graph& a, const Graph& b) {
		return a.n > b.n;
	};
	sort(atoms.begin(), atoms.end(), cmp);
	pair<int, vector<pair<int, int> > > sol = {lowerBound, {}};
	for (int i = 0; i < (int)atoms.size(); i++) {
		pair<int, vector<pair<int, int> > > tSol;
		if (i+1 < (int)atoms.size()) {
			tSol = SolveAtom::solveAtom(atoms[i], 1, sol.F, min(upperBound, fillWidth-1), atoms[i].mapTableInto(vertMap));
		}
		else {
			tSol = SolveAtom::solveAtom(atoms[i], 1, sol.F, min(upperBound, fillWidth-1), atoms[i].mapTableInto(vertMap));
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

pair<ld, vector<pair<int, int> > > solveFracGraph(const Graph& G, ld lowerBound, ld okWidth, ld upperBound, const vector<int>& vertMap) {
	assert(lowerBound > -eps);
	assert(okWidth < upperBound + eps);
	if (lowerBound > upperBound) return {-1, FAILSOLUTION};
	okWidth = max(okWidth, lowerBound);
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	ld fillWidth = MCS::getFracHyperTreeWidth(gFill, vertMap);
	assert(fillWidth > lowerBound - eps);
	
	if (fillWidth < okWidth + eps) return {fillWidth, fill.F};
	
	vector<Graph> atoms = MCS::getAtoms(G, gFill, fill).S;
	
	auto cmp = [](const Graph& a, const Graph& b) {
		return a.n > b.n;
	};
	sort(atoms.begin(), atoms.end(), cmp);
	pair<ld, vector<pair<int, int> > > sol = {lowerBound, {}};
	for (int i = 0; i < (int)atoms.size(); i++) {
		pair<ld, vector<pair<int, int> > > tSol;
		if (i+1 < (int)atoms.size()) {
			tSol = SolveAtom::solveFracAtom(atoms[i], 1, sol.F, min(upperBound, fillWidth), atoms[i].mapTableInto(vertMap));
		}
		else {
			tSol = SolveAtom::solveFracAtom(atoms[i], 1, sol.F, min(upperBound, fillWidth), atoms[i].mapTableInto(vertMap));
		}
		if (tSol.S == FAILSOLUTION) {
			if (fillWidth < upperBound + eps) {
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
	assert(sol.F < upperBound + eps);
	return sol;
}

}

class HyperIO {
	vector<string> vertS;
	
public:
	vector<vector<int> > readInput(istream& in) {
		map<string, int> vMap;
		vector<vector<int> > hEdges;
		string inS;
		while (getline(in, inS)) {
			for (char& c : inS) {
				if (c == ',') c = ' ';
			}
			stringstream ss;
			ss<<inS;
			string ve;
			vector<int> te;
			while (ss>>ve) {
				if (!vMap.count(ve)) {
					vMap[ve] = vertS.size();
					vertS.push_back(ve);
				}
				te.push_back(vMap[ve]);
			}
			sort(te.begin(), te.end());
			te.erase(unique(te.begin(), te.end()), te.end());
			if (te.size() < 2) continue;
			hEdges.push_back(te);
		}
		return hEdges;
	}
	
	void printFill(ostream& out, vector<pair<int, int> > edges) {
		for (auto e : edges) {
			assert(e.F>=0&&e.F<(int)vertS.size()&&e.S>=0&&e.S<(int)vertS.size());
			out<<vertS[e.F]<<" "<<vertS[e.S]<<endl;
		}
	}
};

void solve() {
	Timer timer;
	timer.start();
	HyperIO io;
	vector<vector<int> > hEdges = io.readInput(cin);
	set<pair<int, int> > pEdges;
	for (auto& he : hEdges) {
		for (int i = 0; i < (int)he.size(); i++) {
			for (int ii = i+1; ii < (int)he.size(); ii++) {
				pEdges.insert({he[i], he[ii]});
			}
		}
	}
	Graph pGraph(pEdges);
	vector<vector<int> > hes = hEdges;
	for (auto& he : hes) {
		for (int& x : he) {
			x = pGraph.mapInto(x);
		}
		sort(he.begin(), he.end());
	}
	SetCover::init(hes);
	vector<int> vert(pGraph.n);
	for (int i = 0; i < pGraph.n; i++) {
		vert[i] = i;
	}
	auto sol = SolveGraph::solveGraph(pGraph, 0, 0, (int)hes.size(), vert);
	sol.S = pGraph.mapBack(sol.S);
	io.printFill(cout, sol.S);
	cerr<<"Generalized hypertreewidth: "<<sol.F<<endl;
	cerr<<"*Fill edges: "<<sol.S.size()<<endl;
	cerr<<"Time "<<timer.getTime().count()<<endl;
}

void solveFrac() {
	Timer timer;
	timer.start();
	HyperIO io;
	vector<vector<int> > hEdges = io.readInput(cin);
	set<pair<int, int> > pEdges;
	for (auto& he : hEdges) {
		for (int i = 0; i < (int)he.size(); i++) {
			for (int ii = i+1; ii < (int)he.size(); ii++) {
				pEdges.insert({he[i], he[ii]});
			}
		}
	}
	Graph pGraph(pEdges);
	vector<vector<int> > hes = hEdges;
	for (auto& he : hes) {
		for (int& x : he) {
			x = pGraph.mapInto(x);
		}
		sort(he.begin(), he.end());
	}
	SetCover::init(hes);
	vector<int> vert(pGraph.n);
	for (int i = 0; i < pGraph.n; i++) {
		vert[i] = i;
	}
	auto sol = SolveGraph::solveFracGraph(pGraph, 0, 0, (ld)hes.size(), vert);
	sol.S = pGraph.mapBack(sol.S);
	io.printFill(cout, sol.S);
	cerr<<setprecision(8)<<"Fractional hypertreewidth: "<<sol.F<<endl;
	cerr<<"*Fill edges: "<<sol.S.size()<<endl;
	cerr<<"Time "<<timer.getTime().count()<<endl;
}
}