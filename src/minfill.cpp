#include "MCS.hpp"
#include "minfill.hpp"
#include "graph.hpp"
#include "PMC.hpp"
#include "global.hpp"

#include <vector>
#include <cassert>
#include <algorithm>
#include <queue>
#include <iostream>

#define F first
#define S second

using namespace std;

namespace MinFill {

namespace BackTrack {

vector<pair<int, int> > searchNaive(Graph& G, int lowerBound, int upperBound, int step, int depth) {
	assert(lowerBound >= 0);
	if (upperBound < lowerBound) return FAILSOLUTION;
	vector<int> cycle = MCS::getChordlessCycle(G);
	if (cycle.size() == 0) {
		return {};
	}
	lowerBound = max(lowerBound, (int)cycle.size()-3);
	if (lowerBound > upperBound) return FAILSOLUTION;
	
	if (step == 3) {
		auto sol = SolveGraph::solveGraph(G, lowerBound, upperBound, depth);
		return sol;
	}
	
	G.addEdge(cycle[0], cycle[2]);
	auto sol1 = searchNaive(G, max(0, lowerBound-1), upperBound-1, step+1, depth+1);
	G.popEdge(cycle[0], cycle[2]);
	if (sol1 != FAILSOLUTION) {
		sol1.push_back({cycle[0], cycle[2]});
		assert((int)sol1.size() <= upperBound && (int)sol1.size() >= lowerBound);
		upperBound = min(upperBound, (int)sol1.size()-1);
		if ((int)sol1.size() == lowerBound) return sol1;
	}
	
	for (int i = 3; i < (int)cycle.size(); i++) {
		G.addEdge(cycle[1], cycle[i]);
		auto sol0 = searchNaive(G, max(0, lowerBound-1), upperBound-1, step+1, depth+1);
		G.popEdge(cycle[1], cycle[i]);
		
		if (sol0 != FAILSOLUTION) {
			sol0.push_back({cycle[1], cycle[i]});
			assert((int)sol0.size() <= upperBound && (int)sol0.size() >= lowerBound);
			upperBound = min(upperBound, (int)sol0.size()-1);
			if ((int)sol0.size() == lowerBound) return sol0;
			if (sol1 != FAILSOLUTION) {
				assert(sol0.size() < sol1.size());
			}
			sol1 = sol0;
		}
	}
	return sol1;
}
vector<pair<int, int> > solve(Graph G, int lowerBound, int upperBound, int depth) {
	assert(lowerBound >= 0);
	if (upperBound < lowerBound) return FAILSOLUTION;
	auto sol = searchNaive(G, lowerBound, upperBound, 0, depth);
	return sol;
}

}

namespace SolveAtom {

int del1[maxN];
int d[maxN];
int from[maxN];
int nb1[maxN];
int bfs[maxN];
int tmp[maxN];
vector<int> findCycle2(const Graph& G, int x, int y, vector<int>& removedFromB) {
	int bfsBack = 0;
	int tmpBack = 0;
	bfs[bfsBack++] = y;
	d[y] = 1;
	int i2 = 1;
	for (int ny : G.g[y]) {
		if (!del1[ny] && !removedFromB[ny]) {
			d[ny] = 1;
			from[ny] = y;
			bfs[bfsBack++] = ny;
		}
		else if (ny != x) {
			d[ny] = 1;
			tmp[tmpBack++] = ny;
		}
	}
	bool fo = false;
	while (i2 < bfsBack) {
		int t = bfs[i2++];
		for (int nt : G.g[t]) {
			if (d[nt]) continue;
			if (del1[nt]) continue;
			if (nt == x) {
				fo = true;
				from[x] = t;
				bfs[bfsBack++] = nt;
				break;
			}
			d[nt] = d[t]+1;
			from[nt] = t;
			bfs[bfsBack++] = nt;
		}
		if (fo) break;
	}
	vector<int> ret;
	if (fo) {
		while (x != y) {
			ret.push_back(x);
			x = from[x];
		}
		ret.push_back(y);
	}
	for (int i = 0; i < bfsBack; i++) {
		d[bfs[i]] = 0;
		from[bfs[i]] = 0;
	}
	for (int i = 0; i < tmpBack; i++) {
		d[tmp[i]] = 0;
		from[tmp[i]] = 0;
	}
	return ret;
}
int findPath(const Graph& G, int x, int y) {
	int bfsBack = 0;
	bfs[bfsBack++] = x;
	d[x] = 1;
	int i2 = 0;
	bool fo = false;
	while (i2 < bfsBack) {
		int t = bfs[i2++];
		for (int nt : G.g[t]) {
			if (d[nt]) continue;
			if (del1[nt]) continue;
			if (nt == y) {
				fo = true;
				break;
			}
			d[nt] = 1;
			bfs[bfsBack++] = nt;
		}
		if (fo) break;
	}
	for (int i = 0; i < bfsBack; i++) {
		d[bfs[i]] = 0;
	}
	return fo;
}
// Find the best solution that is at most upperBound. There is no solution better than lowerBound.
vector<pair<int, int> > solveAtom1(const Graph& G, int lowerBound, int upperBound, int depth) {
	assert(lowerBound >= 0);
	if (lowerBound > upperBound) return FAILSOLUTION;
	auto fill = MCS::MCS_MP(G);
	if (fill.F.size() <= 1 || (int)fill.F.size() <= lowerBound) {
		assert((int)fill.F.size() >= lowerBound);
		if ((int)fill.F.size() > upperBound) return FAILSOLUTION;
		return fill.F;
	}
	lowerBound = max(lowerBound, 1);
	upperBound = min(upperBound, (int)fill.F.size());
	
	// Kernelization part 1, extracting independent chordless cycles
	vector<int> removedFromB(G.n);
	vector<int> A;
	int cc = 0;
	int startVertex = 0;
	while (startVertex < G.n) {
		vector<int> cycle = MCS::getChordlessCycle(G, startVertex, removedFromB);
		if (cycle.size() == 0) {
			break;
		}
		cc += (int)cycle.size()-3;
		for (int x : cycle) {
			assert(removedFromB[x] == 0);
			removedFromB[x] = 1;
			A.push_back(x);
		}
		while (startVertex < G.n && removedFromB[startVertex]) startVertex++;
	}
	lowerBound = max(lowerBound, cc);
	
	if (lowerBound > upperBound) return FAILSOLUTION;
	if (lowerBound >= (int)fill.F.size()) {
		assert((int)fill.F.size() >= lowerBound);
		return fill.F;
	}
	// Kernelization part 2, extracting related chordless cycles
	for (int it = 0; it < (int)A.size(); it++) {
		int x = A[it];
		for (int y : G.g[x]) {
			if (removedFromB[y]) continue;
			for (int nx : G.g[x]) {
				del1[nx]++;
			}
			for (int ny : G.g[y]) {
				if (del1[ny]) del1[ny]++;
			}
			for (int nx : G.g[x]) {
				if (del1[nx] == 2) {
					del1[nx] = 1;
				}
				else {
					del1[nx] = 0;
				}
			}
			vector<int> cycle = findCycle2(G, x, y, removedFromB);
			if (cycle.size() > 0) {
				lowerBound = max(lowerBound, (int)cycle.size()-3);
				vector<vector<int> > mPaths;
				bool bf = false;
				for (int i = 0; i < (int)cycle.size(); i++) {
					if (removedFromB[cycle[i]] == 0) {
						mPaths.push_back(vector<int>());
						while (i < (int)cycle.size() && removedFromB[cycle[i]] == 0) {
							mPaths.back().push_back(cycle[i]);
							i++;
						}
						i--;
					}
					else {
						bf = true;
					}
				}
				assert(bf);
				if (mPaths.size() == 1) {
					assert(mPaths[0].size() <= cycle.size()-1);
					if (mPaths[0].size() == cycle.size()-1) {
						cc += (int)cycle.size()-3;
					}
					else {
						cc += (int)mPaths[0].size()-1;
					}
				}
				else {
					int su = 0;
					int ma = 0;
					for (auto& mPath : mPaths) {
						su += (int)mPath.size()-1;
						ma = max(ma, (int)mPath.size()-1);
					}
					cc += max(su/2, ma);
				}
				for (auto& mPath : mPaths) {
					for (int v : mPath) {
						assert(removedFromB[v] == 0);
						removedFromB[v] = 1;
						A.push_back(v);
					}
				}
				lowerBound = max(lowerBound, cc);
			}
			for (int nx : G.g[x]) {
				del1[nx] = 0;
			}
			for (int ny : G.g[y]) {
				del1[ny] = 0;
			}
		}
	}
	
	if (lowerBound > upperBound) return FAILSOLUTION;
	if (lowerBound >= (int)fill.F.size()) {
		assert(lowerBound <= (int)fill.F.size());
		return fill.F;
	}
	// Kernelization part 3, finding essential edges
	vector<pair<pair<int, int>, vector<int> > > As;
	for (int it1 = 0; it1 < (int)A.size(); it1++) {
		for (int it2 = it1+1; it2 < (int)A.size(); it2++) {
			int y = A[it1];
			int z = A[it2];
			bool f = false;
			for (int ny : G.g[y]) {
				if (ny == z) {
					f = true;
					break;
				}
			}
			if (f) continue;
			for (int ny : G.g[y]) {
				nb1[ny] = 1;
			}
			vector<int> xx;
			for (int nz : G.g[z]) {
				if (nb1[nz]) xx.push_back(nz);
			}
			for (int ny : G.g[y]) {
				nb1[ny] = 0;
			}
			vector<int> AA;
			for (int x : xx) {
				for (int nx : G.g[x]) {
					if (nx != y && nx != z) {
						del1[nx] = 1;
					}
				}
				del1[x] = 1;
				if (findPath(G, y, z)) {
					AA.push_back(x);
				}
				for (int nx : G.g[x]) {
					del1[nx] = 0;
				}
				del1[x] = 0;
			}
			As.push_back({{y, z}, AA});
		}
	}
	auto cmp = [](const pair<pair<int, int>, vector<int> >& a, const pair<pair<int, int>, vector<int> >& b) {
		return a.S.size() < b.S.size();
	};
	sort(As.begin(), As.end(), cmp);
	
	int pIt = 0;
	while (pIt < (int)As.size() && (int)As[pIt].S.size() <= 2*lowerBound) {
		for (int x : As[pIt].S) {
			if (removedFromB[x] == 0) {
				removedFromB[x] = 1;
				A.push_back(x);
			}
		}
		pIt++;
	}
	As.erase(As.begin(), As.begin() + pIt);
	
	if (lowerBound > upperBound) return FAILSOLUTION;
	if (lowerBound >= (int)fill.F.size()) {
		assert(lowerBound <= (int)fill.F.size());
		return fill.F;
	}
	
	int kk = lowerBound;
	if ((int)As.size() > 0) assert((int)As[0].S.size() > 2*kk);
	
	vector<pair<int, int> > kernelIv;
	
	for (int i = 0; i < (int)As.size(); i++) {
		if ((int)As[i].S.size() > 2*kk) {
			kernelIv.push_back({kk, ((int)As[i].S.size()-1)/2});
			assert(((int)As[i].S.size()-1)/2 >= kk);
			kk = ((int)As[i].S.size()-1)/2+1;
		}
	}
	
	if (kk <= upperBound) kernelIv.push_back({kk, upperBound});
	assert(!kernelIv.empty());
	while (kernelIv.back().F > upperBound) {
		kernelIv.pop_back();
		assert(!kernelIv.empty());
	}
	if (kernelIv.back().S > upperBound) {
		kernelIv.back().S = upperBound;
	}
	for (int i = 0; i < (int)kernelIv.size(); i++) {
		assert(kernelIv[i].F <= kernelIv[i].S);
		if (i > 0) {
			assert(kernelIv[i].F == kernelIv[i-1].S+1);
		}
	}
	
	vector<pair<int, int> > solOpt = FAILSOLUTION;
	for (auto iv : kernelIv) {
		assert(iv.F >= lowerBound && iv.S <= upperBound);
		vector<int> AA = A;
		vector<int> inA = removedFromB;
		vector<pair<int, int> > kes;
		vector<pair<int, int> > sol0;
		for (auto& t : As) {
			if ((int)t.S.size() <= 2*iv.S) {
				assert((int)t.S.size() <= 2*iv.F);
				for (int x : t.S) {
					if (!inA[x]) {
						inA[x] = 1;
						AA.push_back(x);
					}
				}
			}
			else {
				kes.push_back(t.F);
				sol0.push_back(t.F);
			}
		}
		for (int x : AA) {
			for (int nx : G.g[x]) {
				if (nx < x && inA[nx]) {
					kes.push_back({nx, x});
				}
			}
		}
		iv.F = max(iv.F, (int)sol0.size());
		if (iv.F > iv.S) {
			lowerBound = iv.S+1;
			continue;
		}
		Graph kernel(kes);
		assert(iv.F >= (int)sol0.size());
		assert(kernel.n <= G.n);
		bool hd12 = false;
		for (int i = 0; i < kernel.n; i++) {
			if (kernel.g[i].size() == 1 || kernel.g[i].size() == 2) {
				hd12 = true;
			}
		}
		vector<pair<int, int> > sol;
		if (sol0.size() > 0 || kernel.n < G.n || kernel.m < G.m) {
			sol = SolveGraph::solveGraph(kernel, iv.F-(int)sol0.size(), iv.S-(int)sol0.size(), depth);
		}
		else {
			assert(hd12 == false);
			sol = PMC::solve(G, (int64_t)(iv.F-(int)sol0.size()), (int64_t)(iv.S-(int)sol0.size()), PROBLEM_MINFILL);
		}
		if (sol != FAILSOLUTION) {
			sol = kernel.mapBack(sol);
			sol.insert(sol.end(), sol0.begin(), sol0.end());
			assert((int)sol.size() >= lowerBound && (int)sol.size() <= iv.S);
			return sol;
		}
		else {
			lowerBound = iv.S+1;
		}
	}
	return FAILSOLUTION;
}
vector<pair<int, int> > greedyDegree2(Graph& G) {
	queue<int> d1;
	queue<int> d2;
	for (int i = 0; i < G.n; i++) {
		assert(G.g[i].size() > 1);
		if (G.g[i].size() == 2) {
			d2.push(i);
		}
	}
	vector<pair<int, int> > sol;
	while (!d1.empty() || !d2.empty()) {
		while (!d1.empty()) {
			int x = d1.front();
			d1.pop();
			assert(G.g[x].size() <= 1);
			if (G.g[x].size() == 0) continue;
			int n1 = G.g[x][0];
			G.remEdge(x, n1);
			if (G.g[n1].size() == 1) d1.push(n1);
			if (G.g[n1].size() == 2) d2.push(n1);
		}
		if (!d2.empty()) {
			int x = d2.front();
			d2.pop();
			assert(G.g[x].size() == 0 || G.g[x].size() == 2);
			if (G.g[x].size() == 0) continue;
			int n1 = G.g[x][0];
			int n2 = G.g[x][1];
			if (!G.hasEdge(n1, n2)) {
				sol.push_back({n1, n2});
				G.addEdge(n1, n2);
			}
			G.remEdge(x, n1);
			G.remEdge(x, n2);
			if (G.g[n1].size() == 1) d1.push(n1);
			if (G.g[n1].size() == 2) d2.push(n1);
			if (G.g[n2].size() == 1) d1.push(n2);
			if (G.g[n2].size() == 2) d2.push(n2);
		}
	}
	return sol;
}
pair<int, pair<int, int> > nbClique(Graph& G) {
	int minF = (int)1e9;
	for (int x = 0; x < G.n; x++) {
		assert(G.g[x].size() > 2);
		pair<int, int> fo = {-1, -1};
		bool fail = false;
		bool brk = false;
		int cnt = 0;
		for (int i = 0; i < (int)G.g[x].size(); i++) {
			for (int ii = i+1; ii < (int)G.g[x].size(); ii++) {
				if (!G.hasEdge(G.g[x][i], G.g[x][ii])) {
					cnt++;
					if (fo.F != -1) {
						fail = true;
						if (cnt >= minF) {
							brk = true;
							break;
						}
					}
					else {
						fo = {G.g[x][i], G.g[x][ii]};
					}
				}
			}
			if (brk) break;
		}
		if (!fail) {
			assert(fo.F != -1);
			return {-1, fo};
		}
		minF = min(minF, cnt);
	}
	return {minF, {-1, -1}};
}
// Find the best solution that is at most upperBound. There is no solution better than lowerBound.
vector<pair<int, int> > solveAtom0(Graph& G, int lowerBound, int upperBound, int depth) {
	assert(lowerBound >= 0);
	if (lowerBound > upperBound) return FAILSOLUTION;
	if (G.n <= 3) return {};
	if ((int64_t)G.n*(int64_t)(G.n-1) == (int64_t)G.m*2) return {};
	if (upperBound <= 0) {
		auto cycle = MCS::getChordlessCycle(G);
		if (cycle.size() == 0) return {};
		else return FAILSOLUTION;
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
		auto sol = SolveGraph::solveGraph(G, max(0, lowerBound - (int)sol0.size()), upperBound - (int)sol0.size(), depth);
		if (sol == FAILSOLUTION) return FAILSOLUTION;
		sol.insert(sol.end(), sol0.begin(), sol0.end());
		assert((int)sol.size() <= upperBound);
		return sol;
	}
	auto nbC = nbClique(G);
	if (nbC.S.F != -1) {
		auto sol0 = nbC.S;
		G.addEdge(sol0.F, sol0.S);
		auto sol = SolveGraph::solveGraph(G, max(0, lowerBound-1), upperBound-1, depth);
		if (sol == FAILSOLUTION) return FAILSOLUTION;
		sol.push_back(sol0);
		assert((int)sol.size() <= upperBound);
		return sol;
	}
	else {
		assert(nbC.F >= 2);
		lowerBound = max(lowerBound, nbC.F);
		if (lowerBound > upperBound) return FAILSOLUTION;
	}
	auto ret = solveAtom1(G, lowerBound, upperBound, depth);
	return ret;
}

}
namespace SolveGraph {

// Find the best solution that is at most upperBound. There is no solution better than lowerBound.
vector<pair<int, int> > solveGraph(const Graph& G, int lowerBound, int upperBound, int depth) {
	assert(lowerBound >= 0);
	if (lowerBound > upperBound) return FAILSOLUTION;
	auto fill = MCS::MCS_MP(G);
	if (fill.F.size() <= 1 || (int)fill.F.size() <= lowerBound) {
		assert((int)fill.F.size() >= lowerBound);
		if ((int)fill.F.size() > upperBound) return FAILSOLUTION;
		return fill.F;
	}
	lowerBound = max(lowerBound, 1);
	upperBound = min(upperBound, (int)fill.F.size());
	if (lowerBound > upperBound) return FAILSOLUTION;
	
	Graph gFill = G;
	gFill.addEdges(fill.F);
	vector<Graph> atoms = MCS::getAtoms(G, gFill, fill).S;
	auto cmp = [](const Graph& a, const Graph& b) {
		return a.n < b.n;
	};
	sort(atoms.begin(), atoms.end(), cmp);
	vector<pair<int, int> > sol;
	for (int i = 0; i < (int)atoms.size(); i++) {
		vector<pair<int, int> > tSol;
		if (i+1 < (int)atoms.size()) {
			tSol = SolveAtom::solveAtom0(atoms[i], 0, upperBound - (int)sol.size(), depth);
		}
		else {
			tSol = SolveAtom::solveAtom0(atoms[i], max(0, lowerBound - (int)sol.size()), upperBound - (int)sol.size(), depth);
		}
		if (tSol == FAILSOLUTION) return FAILSOLUTION;
		tSol = atoms[i].mapBack(tSol);
		sol.insert(sol.end(), tSol.begin(), tSol.end());
	}
	assert((int)sol.size() >= lowerBound);
	assert((int)sol.size() <= upperBound);
	return sol;
}

}

}