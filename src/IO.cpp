#include "IO.hpp"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>

#define F first
#define S second

using namespace std;

int IO::nTokens(string s) const {
	stringstream ss;
	ss<<s;
	int r=0;
	while (ss>>s) {
		r++;
	}
	return r;
}

int IO::readInput(istream& in, set<pair<int, int> >& edges) {
	ios_base::sync_with_stdio(0);
	string inputLine;
	bool dimacs = false;
	paceTreeWidth = false;
	paceVertices = 0;
	assert(edges.size() == 0);
	while (getline(in, inputLine)) {
		if (inputLine[0] == '#') continue;
		if (inputLine[0] == 'c' && (dimacs || paceTreeWidth) ) continue;
		if (inputLine.size() >= 7 && inputLine.substr(0, 7) == "p edge " && nTokens(inputLine) >= 3) {
			dimacs = true;
			edges.clear();
			continue;
		}
		if (inputLine.size() >= 6 && inputLine.substr(0, 6) == "p col " && nTokens(inputLine) >= 3) {
			dimacs = true;
			edges.clear();
			continue;
		}
		if (inputLine.size() >= 5 && inputLine.substr(0, 5) == "p tw " && nTokens(inputLine) >= 3) {
			edges.clear();
			stringstream ss;
			ss<<inputLine;
			string temp;
			ss>>temp;
			ss>>temp;
			ss>>paceVertices;
			paceTreeWidth = true;
			continue;
		}
		if (inputLine.size() >= 2 && inputLine.substr(0, 2) == "p " && nTokens(inputLine) >= 3) {
			edges.clear();
			dimacs = true;
			continue;
		}
		string vA,vB;
		stringstream ss;
		ss<<inputLine;
		if (inputLine[0] == 'e' && dimacs && nTokens(inputLine) == 3) {
			char e;
			if (!(ss>>e>>vA>>vB)) {
				cerr<<"Unexpected input "<<inputLine<<endl;
				edges.clear();
				continue;
			}
		}
		else if (nTokens(inputLine) == 2 && !dimacs) {
			if (!(ss>>vA>>vB)) {
				cerr<<"Unexpected input "<<inputLine<<endl;
				edges.clear();
				continue;
			}
		}
		else {
			continue;
		}
		if (vA == vB) continue;
		if (vertices.count(vA) == 0) {
			verticeLabels.push_back(vA);
			vertices[vA] = verticeLabels.size();
		}
		if (vertices.count(vB) == 0) {
			verticeLabels.push_back(vB);
			vertices[vB] = verticeLabels.size();
		}
		int a = vertices[vA];
		int b = vertices[vB];
		if (a > b) swap(a, b);
		edges.insert({a, b});
	}
	return vertices.size();
}
string IO::vertexId(int id) const {
	return verticeLabels[id-1];
}
void IO::printFill(ostream& out, const vector<pair<int, int> >& edges) const {
	for (auto e : edges) {
		out<<vertexId(e.F)<<" "<<vertexId(e.S)<<'\n';
	}
	out.flush();
}
void IO::printTreeDecomposition(ostream& out, pair<vector<vector<int> >, vector<pair<int, int> > > treeDecomposition) {
	if (paceTreeWidth) {
		set<string> hasVertices;
		for (auto& bag : treeDecomposition.F) {
			for (int x : bag) {
				hasVertices.insert(vertexId(x));
			}
		}
		for (int i = 1; i <= paceVertices; i++) {
			if (!hasVertices.count(to_string(i))) {
				verticeLabels.push_back(to_string(i));
				vertices[to_string(i)] = verticeLabels.size();
				treeDecomposition.F.push_back({(int)verticeLabels.size()});
				
				if (treeDecomposition.F.size() > 1) {
					treeDecomposition.S.push_back({0, (int)treeDecomposition.F.size()-1});
				}
			}
		}
	}
	int maxBag = 0;
	for (auto& bag : treeDecomposition.F) {
		maxBag = max(maxBag, (int)bag.size());
	}
	if (paceTreeWidth) {
		out<<"s td "<<treeDecomposition.F.size()<<" "<<maxBag<<" "<<paceVertices<<'\n';
	}
	else {
		out<<"s td "<<treeDecomposition.F.size()<<" "<<maxBag<<" "<<verticeLabels.size()<<'\n';
	}
	for (int i = 0; i < (int)treeDecomposition.F.size(); i++) {
		out<<"b "<<i+1;
		for (int x : treeDecomposition.F[i]) {
			out<<" "<<vertexId(x);
		}
		out<<'\n';
	}
	for (auto e : treeDecomposition.S) {
		out<<e.F+1<<" "<<e.S+1<<'\n';
	}
	out.flush();
}