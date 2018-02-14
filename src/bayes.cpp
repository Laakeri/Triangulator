#include "bayes.hpp"
#include "global.hpp"
#include "graph.hpp"
#include "MCS.hpp"
#include "PMC.hpp"

#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <cassert>
#include <set>
#include <algorithm>

#define F first
#define S second

using namespace std;

namespace Bayes {

namespace SolveAtom {
pair<int64_t, vector<pair<int, int> > > solveAtom(const Graph& G, int64_t lowerBound, int64_t upperBound, vector<int64_t> tableSizes) {
	auto fill = MCS::MCS_MP(G);
	
	if (fill.F.size() <= 0) {
		Graph gFill = G;
		gFill.addEdges(fill.F);
		int64_t fillWidth = MCS::getBayesWidth(gFill, tableSizes);
		assert(fillWidth >= lowerBound);
		if (fillWidth > upperBound) return {-1, FAILSOLUTION};
		return {fillWidth, fill.F};
	}
	
	PMC::initTableSz(tableSizes);
	auto fillO = PMC::solve(G, lowerBound, upperBound, PROBLEM_BAYES);
	if (fillO == FAILSOLUTION) return {-1, FAILSOLUTION};
	Graph gFillO = G;
	gFillO.addEdges(fillO);
	int width = MCS::getBayesWidth(gFillO, tableSizes);
	assert(width >= lowerBound);
	assert(width <= upperBound);
	return {width, fillO};
}
}

namespace SolveGraph {

// Find a fill-in with width at most okWidth or the best solution that is at most upperBound. There is no solution better than lowerBound.
pair<int64_t, vector<pair<int, int> > > solveGraph(const Graph& G, int64_t lowerBound, int64_t upperBound, vector<int64_t> tableSizes) {
	assert(lowerBound >= 0);
	assert(lowerBound <= upperBound);
	assert((int)tableSizes.size() == G.n);
	for (int i = 0; i < G.n; i++) {
		for (int nx : G.g[i]) {
			lowerBound = max(lowerBound, tableSizes[i]*tableSizes[nx]);
		}
	}
	
	auto fill = MCS::MCS_MP(G);
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	
	int64_t fillWidth = MCS::getBayesWidth(gFill, tableSizes);
	
	assert(fillWidth >= lowerBound);
	if (fillWidth <= lowerBound) return {fillWidth, fill.F};
	
	if (fill.F.size() <= 0) {
		if (fillWidth > upperBound) return {-1, FAILSOLUTION};
		return {fillWidth, fill.F};
	}
	
	vector<Graph> atoms = MCS::getAtoms(G, gFill, fill).S;
	auto cmp = [](const Graph& a, const Graph& b) {
		return a.n > b.n;
	};
	sort(atoms.begin(), atoms.end(), cmp);
	pair<int64_t, vector<pair<int, int> > > sol = {0, {}};
	for (int i = 0; i < (int)atoms.size(); i++) {
		pair<int64_t, vector<pair<int, int> > > tSol;
		if (i+1 < (int)atoms.size()) {
			tSol = SolveAtom::solveAtom(atoms[i], 0, upperBound - sol.F, atoms[i].mapTableInto(tableSizes));
		}
		else {
			tSol = SolveAtom::solveAtom(atoms[i], max((int64_t)0, lowerBound - sol.F), upperBound - sol.F, atoms[i].mapTableInto(tableSizes));
		}
		if (tSol.S == FAILSOLUTION) {
			if (fillWidth <= upperBound) {
				return {fillWidth, fill.F};
			}
			else {
				return {-1, FAILSOLUTION};
			}
		}
		sol.F += tSol.F;
		if (sol.F > upperBound) return {-1, FAILSOLUTION};
		tSol.S = atoms[i].mapBack(tSol.S);
		sol.S.insert(sol.S.end(), tSol.S.begin(), tSol.S.end());
	}
	assert(sol.F <= upperBound);
	return sol;
}
}


class BayesIO {
public:
	set<pair<int, int> > dirEdges;
	map<string, int> variables;
	vector<int64_t> tableSizes;
	int vars = 0;
	
	void readInput(istream& in) {
		string inS;
		stringstream ss;
		int qde=0;
		while (getline(in, inS)) {
			for (int i = 0; i < (int)inS.size(); i++) {
				if (inS[i] == '/' && i+1 < (int)inS.size() && inS[i+1] == '/') break;
				if (inS[i] == '"') qde = 1-qde;
				
				if (qde == 0) {
					if (inS[i] == '[') ss<<" ";
					if (inS[i] == ']') ss<<" ";
					if (inS[i] == '}') ss<<" ";
					if (inS[i] == '{') ss<<" ";
					if (inS[i] == ')') ss<<" ";
					if (inS[i] == '(') ss<<" ";
					if (inS[i] == ',') inS[i] = ' ';
				}
				ss<<inS[i];
				if (qde == 0) {
					if (inS[i] == '[') ss<<" ";
					if (inS[i] == ']') ss<<" ";
					if (inS[i] == '}') ss<<" ";
					if (inS[i] == '{') ss<<" ";
					if (inS[i] == ')') ss<<" ";
					if (inS[i] == '(') ss<<" ";
				}
			}
			ss<<" ";
		}
		int de = 0;
		int tvar = -1;
		int format = 0;
		while (ss>>inS) {
			if (inS == "{") de++;
			else if (inS == "}") {
				de--;
				if (de == 0) tvar = -1;
			}
			if (format == 0) {
				// .bif
				if (de == 0 && inS == "variable") format = 1;
				// .dsc
				if (de == 0 && inS == "node") format = 2;
				// .net
				if (de == 0 && inS == "net") format = 3;
			}
			if (format == 1) {
				if (de == 0 && inS == "variable") {
					string var;
					ss>>var;
					assert(variables.count(var) == 0);
					tvar = vars++;
					variables[var] = tvar;
				}
				if (de == 0 && inS == "probability") {
					string te;
					ss>>te;
					assert(te == "(");
					string var1;
					ss>>var1;
					ss>>te;
					assert(variables.count(var1));
					if (te == "|") {
						string var;
						while (ss>>var) {
							if (var == ")") break;
							assert(variables.count(var));
							dirEdges.insert({variables[var], variables[var1]});
						}
					}
				}
				if (de == 1 && tvar >= 0 && inS == "type") {
					string type;
					ss>>type;
					assert(type == "discrete");
					ss>>type;
					assert(type == "[");
					int dim;
					ss>>dim;
					assert(dim>0);
					ss>>type;
					assert(type == "]");
					tableSizes.push_back(dim);
					tvar = -1;
				}
			}
			if (format == 2) {
				if (de == 0 && inS == "node") {
					string var;
					ss>>var;
					assert(variables.count(var) == 0);
					tvar = vars++;
					variables[var] = tvar;
				}
				if (de == 0 && inS == "probability") {
					string te;
					ss>>te;
					assert(te == "(");
					string var1;
					ss>>var1;
					ss>>te;
					assert(variables.count(var1));
					if (te == "|") {
						string var;
						while (ss>>var) {
							if (var == ")") break;
							assert(variables.count(var));
							dirEdges.insert({variables[var], variables[var1]});
						}
					}
				}
				if (de == 1 && tvar >= 0 && inS == "type:") {
					string type;
					ss>>type;
					assert(type == "discrete");
					ss>>type;
					assert(type == "[");
					int dim;
					ss>>dim;
					assert(dim>0);
					ss>>type;
					assert(type == "]");
					tableSizes.push_back(dim);
					tvar = -1;
				}
			}
			if (format == 3) {
				if (de == 0 && inS == "node") {
					string var;
					ss>>var;
					assert(variables.count(var) == 0);
					tvar = vars++;
					variables[var] = tvar;
				}
				if (de == 1 && tvar >=0 && inS == "states") {
					string te;
					ss>>te;
					assert(te == "=");
					ss>>te;
					assert(te == "(");
					int dim=0;
					while (1) {
						ss>>te;
						if (te[0] == '"' && te.back() == '"') {
							dim++;
						}
						else {
							assert(te == ")");
							ss>>te;
							assert(te == ";");
							break;
						}
					}
					tableSizes.push_back(dim);
					tvar = -1;
				}
				if (de == 0 && inS == "potential") {
					string te;
					ss>>te;
					assert(te == "(");
					string var1;
					ss>>var1;
					ss>>te;
					assert(variables.count(var1));
					if (te == "|") {
						string var;
						while (ss>>var) {
							if (var == ")") break;
							assert(variables.count(var));
							dirEdges.insert({variables[var], variables[var1]});
						}
					}
				}
			}
		}
	}
	
	void printFill(ostream& out, vector<pair<int, int> > es) {
		vector<string> mp(variables.size());
		for (auto v : variables) {
			assert(v.S>=0&&v.S<(int)variables.size());
			mp[v.S] = v.F;
		}
		for (auto e : es) {
			out<<mp[e.F]<<" "<<mp[e.S]<<endl;
		}
	}
};

void solve() {
	BayesIO io;
	Timer timer;
	timer.start();
	io.readInput(cin);
	assert(io.vars == (int)io.variables.size());
	assert(io.vars == (int)io.tableSizes.size());
	set<pair<int, int> > edges;
	// moralizing
	for (auto e : io.dirEdges) {
		edges.insert(minmax(e.F, e.S));
	}
	for (int i = 0; i < io.vars; i++) {
		vector<int> inc;
		for (int ii = 0; ii < io.vars; ii++) {
			if (i != ii && io.dirEdges.count({ii, i})) {
				inc.push_back(ii);
			}
		}
		for (int a : inc) {
			for (int b : inc) {
				if (a < b) {
					edges.insert({a, b});
				}
			}
		}
	}
	cerr<<"Vars: "<<io.variables.size()<<", Edges: "<<io.dirEdges.size()<<", Moralized edges: "<<edges.size()<<endl;
	Graph G(edges);
	vector<int64_t> tableSizes = G.mapTableInto(io.tableSizes);
	auto sol = SolveGraph::solveGraph(G, 0, 1e18, tableSizes);
	assert(sol.S != FAILSOLUTION);
	G.addEdges(sol.S);
	assert(sol.F == MCS::getBayesWidth(G, tableSizes));
	io.printFill(cout, sol.S);
	cerr<<"Total table size: "<<sol.F<<endl;
	cerr<<"*Fill edges: "<<sol.S.size()<<endl;
	cerr<<"Time "<<timer.getTime().count()<<endl;
}

}