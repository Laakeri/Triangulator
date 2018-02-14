#include "MCS.hpp"
#include "global.hpp"
#include "setcover.hpp"

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

#define F first
#define S second

using namespace std;

typedef long double ld;

namespace MCS {

vector<int> reach[maxN];
vector<int> q[maxN];
int label[maxN];
int rm[maxN];
int rc[maxN];
int ca[maxN];
int lb[maxN];
int tr[maxN];
int d[maxN];
int from[maxN];
int index[maxN];
int f[maxN];
int invOrd[maxN];
int mk[maxN];
// finds a "minimal elimination order": an elimination order that gives a minimal chordalization in O(nm)
pair<vector<pair<int, int> >, vector<pair<int, int> > > MCS_MP(const Graph& G) {
	for (int i = 0; i < G.n; i++) {
		label[i] = 0;
		rm[i] = 0;
		reach[i].clear();
		rc[i] = 0;
	}
	vector<pair<int, int> > fill;
	vector<pair<int, int> > order(G.n);
	int s = -1;
	for (int it = G.n-1; it >= 0; it--) {
		int x = -1;
		for (int i = 0; i < G.n; i++) {
			if (rm[i] == 0 && label[i] >= label[x]) {
				x = i;
			}
		}
		assert(x >= 0 && rm[x] == 0 && label[x] < G.n);
		order[it] = {x, (label[x] <= s)};
		s = label[x];
		for (int i = 0; i < G.n; i++) {
			rc[i] = 0;
		}
		rc[x] = 1;
		rm[x] = 1;
		for (int y : G.g[x]) {
			if (rm[y] == 0) {
				rc[y] = 1;
				reach[label[y]].push_back(y);
			}
		}
		for (int j = 0; j < G.n; j++) {
			while (!reach[j].empty()) {
				int y = reach[j].back();
				reach[j].pop_back();
				for (int z : G.g[y]) {
					if (rm[z] == 0 && rc[z] == 0) {
						rc[z] = 1;
						if (label[z] > j) {
							reach[label[z]].push_back(z);
							label[z]++;
							fill.push_back({x, z});
						}
						else {
							reach[j].push_back(z);
						}
					}
				}
			}
		}
		for (int y : G.g[x]) {
			if (rm[y] == 0) {
				label[y]++;
			}
		}
	}
	return {fill, order};
}
void findCC(int x, const Graph& G, vector<pair<int, int> >& v) {
	if (rc[x]) return;
	rc[x] = 1;
	for (int nx : G.g[x]) {
		if (ca[nx] || !rc[nx]) {
			v.push_back({x, nx});
		}
	}
	for (int nx : G.g[x]) {
		if (!ca[nx] && !rc[nx]) {
			findCC(nx, G, v);
		}
	}
}
// return the atoms of a graph given a minimal elimination order in O(nm), also size of the maximum clique separator
pair<int, vector<Graph> > getAtoms(const Graph& G, const Graph& gFill, const pair<vector<pair<int, int> >, vector<pair<int, int> > >& fill) {
	assert(G.n == gFill.n);
	for (int i = 0; i < G.n; i++) {
		rm[i] = 0;
		rc[i] = 0;
		ca[i] = 0;
	}
	int maxSep = 0;
	vector<Graph> ret;
	for (int i = 0; i < G.n; i++) {
		if (fill.S[i].S) {
			vector<int> cand;
			for (int nx : gFill.g[fill.S[i].F]) {
				if (rm[nx] == 0) {
					cand.push_back(nx);
				}
			}
			for (int x : cand) {
				ca[x] = 1;
			}
			bool no = false;
			for (int x : cand) {
				int fo = 0;
				for (int nx : G.g[x]) {
					if (ca[nx]) fo++;
				}
				if (fo < (int)cand.size()-1) {
					no = true;
				}
			}
			if (!no) {
				maxSep = max(maxSep, (int)cand.size()+1);
				vector<pair<int, int> > es;
				findCC(fill.S[i].F, G, es);
				for (int j = 0; j < (int)cand.size(); j++) {
					for (int jj = j+1; jj < (int)cand.size(); jj++) {
						es.push_back({cand[j], cand[jj]});
					}
				}
				ret.push_back(es);
			}
			for (int x : cand) {
				ca[x] = 0;
			}
		}
		rm[fill.S[i].F] = 1;
	}
	int fo = 0;
	for (int i = 0; i < G.n; i++) {
		if (!rc[i]) {
			vector<pair<int, int> > es;
			findCC(i, G, es);
			ret.push_back(Graph(es));
			fo++;
		}
	}
	assert(fo == 1);
	return {maxSep, ret};
}
// returns an elimination order that is perfect iff the graph is chordal in O(n)
vector<int> MCS_L(const Graph& G, int sv, const vector<int>& rmV) {
	assert(rmV[sv] == 0);
	vector<int> order;
	q[0].push_back(sv);
	int mLb = 0;
	while (mLb >= 0) {
		if (q[mLb].size() == 0) {
			mLb--;
			continue;
		}
		int x = q[mLb].back();
		q[mLb].pop_back();
		if (tr[x]) continue;
		order.push_back(x);
		for (int nx : G.g[x]) {
			if (!rmV[nx] && !tr[nx]) {
				lb[nx]++;
				q[lb[nx]].push_back(nx);
				mLb = max(mLb, lb[nx]);
			}
		}
		tr[x] = 1;
	}
	reverse(order.begin(), order.end());
	for (int x : order) {
		tr[x] = 0;
		lb[x] = 0;
	}
	return order;
}
// returns treewidth of a chordal graph in O(n)
int getTreeWidth(const Graph& G) {
	if (G.n <= 1) return 0;
	vector<int> u(G.n);
	vector<int> order;
	for (int i = 0; i < G.n; i++) {
		if (G.g[i].size() == 0) {
			order.push_back(i);
			u[i] = 1;
		}
		else if (u[i] == 0) {
			auto ord = MCS_L(G, i, u);
			for (int x : ord) {
				assert(u[x] == 0);
				u[x] = 1;
			}
			order.insert(order.end(), ord.begin(), ord.end());
		}
	}
	assert((int)order.size() == G.n);
	for (int i = 0; i < (int)order.size(); i++) {
		invOrd[order[i]] = i;
	}
	int treeWidth = 0;
	for (int i = 0; i < G.n; i++) {
		int x = order[i];
		assert(0<=x&&x<G.n);
		int nb = 0;
		for (int nx : G.g[x]) {
			if (invOrd[nx] > i) {
				nb++;
			}
		}
		treeWidth = max(treeWidth, nb);
	}
	return treeWidth;
}
// returns bayes width of a chordal graph in O(n)
int64_t getBayesWidth(const Graph& G, const vector<int64_t>& tableSizes) {
	if (G.n <= 1) return 0;
	vector<int> u(G.n);
	vector<int> order;
	for (int i = 0; i < G.n; i++) {
		if (G.g[i].size() == 0) {
			order.push_back(i);
			u[i] = 1;
		}
		else if (u[i] == 0) {
			auto ord = MCS_L(G, i, u);
			for (int x : ord) {
				assert(u[x] == 0);
				u[x] = 1;
			}
			order.insert(order.end(), ord.begin(), ord.end());
		}
	}
	assert((int)order.size() == G.n);
	for (int i = 0; i < (int)order.size(); i++) {
		invOrd[order[i]] = i;
	}
	vector<int> fw(G.n);
	for (int i = 0; i < G.n; i++) {
		int x = order[i];
		for (int nx : G.g[x]) {
			if (invOrd[nx] > i) {
				fw[x]++;
			}
		}
	}
	vector<int> noMaxi(G.n);
	int64_t bayesWidth = 0;
	for (int i = 0; i < G.n; i++) {
		int x = order[i];
		assert(0<=x&&x<G.n);
		vector<int> nb;
		for (int nx : G.g[x]) {
			if (invOrd[nx] > i) {
				nb.push_back(nx);
			}
		}
		auto cmp = [&](int a, int b) {
			return invOrd[a] < invOrd[b];
		};
		sort(nb.begin(), nb.end(), cmp);
		for (int ii = 0; ii < (int)nb.size(); ii++) {
			assert(fw[nb[ii]] >= (int)nb.size()-ii-1);
			if (fw[nb[ii]] == (int)nb.size()-ii-1) {
				noMaxi[nb[ii]] = 1;
			}
		}
		if (!noMaxi[x]) {
			int64_t ts = tableSizes[x];
			for (int nx : nb) {
				if ((long double)ts*(long double)tableSizes[nx] > (long double)(maxBWidth)) {
					ts = maxBWidth;
					break;
				}
				ts *= tableSizes[nx];
			}
			bayesWidth += ts;
			if (bayesWidth > maxBWidth) bayesWidth = maxBWidth;
		}
	}
	return bayesWidth;
}
// returns hypertreewidth of a chordal graph in O(n)
int getHyperTreeWidth(const Graph& G, const vector<int>& vertMap) {
	if (G.n <= 1) return 0;
	vector<int> u(G.n);
	vector<int> order;
	for (int i = 0; i < G.n; i++) {
		if (G.g[i].size() == 0) {
			order.push_back(i);
			u[i] = 1;
		}
		else if (u[i] == 0) {
			auto ord = MCS_L(G, i, u);
			for (int x : ord) {
				assert(u[x] == 0);
				u[x] = 1;
			}
			order.insert(order.end(), ord.begin(), ord.end());
		}
	}
	assert((int)order.size() == G.n);
	for (int i = 0; i < (int)order.size(); i++) {
		invOrd[order[i]] = i;
	}
	int treeWidth = 0;
	for (int i = 0; i < G.n; i++) {
		int x = order[i];
		assert(0<=x&&x<G.n);
		vector<int> nb = {vertMap[x]};
		for (int nx : G.g[x]) {
			if (invOrd[nx] > i) {
				nb.push_back(vertMap[nx]);
			}
		}
		auto setCoverSol = SetCover::solve(nb, treeWidth, G.n);
		assert(setCoverSol.size() > 0);
		assert(setCoverSol[0] != -1);
		treeWidth = max(treeWidth, (int)setCoverSol.size());
	}
	return treeWidth;
}
// returns fractional hyper treewidth of a chordal graph in O(n)
ld getFracHyperTreeWidth(const Graph& G, const vector<int>& vertMap) {
	if (G.n <= 1) return 0;
	vector<int> u(G.n);
	vector<int> order;
	for (int i = 0; i < G.n; i++) {
		if (G.g[i].size() == 0) {
			order.push_back(i);
			u[i] = 1;
		}
		else if (u[i] == 0) {
			auto ord = MCS_L(G, i, u);
			for (int x : ord) {
				assert(u[x] == 0);
				u[x] = 1;
			}
			order.insert(order.end(), ord.begin(), ord.end());
		}
	}
	assert((int)order.size() == G.n);
	for (int i = 0; i < (int)order.size(); i++) {
		invOrd[order[i]] = i;
	}
	ld treeWidth = 0;
	for (int i = 0; i < G.n; i++) {
		int x = order[i];
		assert(0<=x&&x<G.n);
		vector<int> nb = {vertMap[x]};
		for (int nx : G.g[x]) {
			if (invOrd[nx] > i) {
				nb.push_back(vertMap[nx]);
			}
		}
		ld setCoverSol = SetCover::solveFrac(nb);
		assert(setCoverSol > 0);
		treeWidth = max(treeWidth, setCoverSol);
	}
	return treeWidth;
}
// returns the tree decomposition of a chordal graph in O(n)
pair<vector<vector<int> >, vector<pair<int, int> > > getTreeDecomposition(const Graph& G) {
	vector<int> u(G.n);
	vector<int> order;
	for (int i = 0; i < G.n; i++) {
		if (G.g[i].size() == 0) {
			order.push_back(i);
			u[i] = 1;
		}
		else if (u[i] == 0) {
			auto ord = MCS_L(G, i, u);
			for (int x : ord) {
				assert(u[x] == 0);
				u[x] = 1;
			}
			order.insert(order.end(), ord.begin(), ord.end());
		}
	}
	assert((int)order.size() == G.n);
	for (int i = 0; i < (int)order.size(); i++) {
		invOrd[order[i]] = i;
	}
	vector<vector<int> > bags(G.n);
	vector<pair<int, int> > edges;
	for (int i = 0; i < G.n; i++) {
		int x = order[i];
		assert(0<=x&&x<G.n);
		bags[i].push_back(x);
		int neB = -1;
		for (int nx : G.g[x]) {
			if (invOrd[nx] > i) {
				bags[i].push_back(nx);
				if (neB == -1 || invOrd[nx] < neB) {
					neB = invOrd[nx];
				}
			}
		}
		if (i < G.n-1) {
			if (neB != -1) {
				edges.push_back({i, neB});
			}
			else {
				edges.push_back({i, i+1});
			}
		}
	}
	assert((int)bags.size() == G.n);
	if (G.n > 0) assert((int)edges.size() == G.n-1);
	return {bags, edges};
}
vector<int> findShortestPath(const Graph& G, int a, int b, const vector<int>& rmV) {
	assert(a != b);
	d[a] = 1;
	vector<int> bfs;
	bfs.push_back(a);
	int i2 = 0;
	bool fo = false;
	while (i2 < (int)bfs.size()) {
		int x = bfs[i2++];
		for (int nx : G.g[x]) {
			if (d[nx]) continue;
			if (rmV[nx]) continue;
			if (x == a && nx == b) continue;
			if (nx == b) {
				from[nx] = x;
				fo = true;
				bfs.push_back(nx);
				break;
			}
			d[nx] = d[x]+1;
			from[nx] = x;
			bfs.push_back(nx);
		}
		if (fo) break;
	}
	vector<int> ret;
	if (fo) {
		while (b != a) {
			ret.push_back(b);
			b = from[b];
		}
		ret.push_back(a);
	}
	for (int x : bfs) {
		d[x] = 0;
		from[x] = 0;
	}
	return ret;
}
// returns a chordless cycle in G in component of sv in O(n) time
vector<int> getChordlessCycle(const Graph& G, int sv, vector<int>& rmV) {
	vector<int> order = MCS_L(G, sv, rmV);
	for (int i = 0; i < (int)order.size(); i++) {
		invOrd[order[i]] = i;
	}
	vector<int> ret;
	for (int i = 0; i < (int)order.size(); i++) {
		int w = order[i];
		f[w] = w;
		index[w] = i;
		for (int v : G.g[w]) {
			if (invOrd[v] >= i || rmV[v]) continue;
			index[v] = i;
			if (f[v] == v) f[v] = w;
		}
		for (int v : G.g[w]) {
			if (invOrd[v] >= i || rmV[v]) continue;
			if (index[f[v]] < i) {
				for (int nx : G.g[v]) {
					if (!rmV[nx]) mk[nx] = 1;
				}
				for (int nx : G.g[w]) {
					if (mk[nx]) rmV[nx] = 1;
				}
				ret = findShortestPath(G, v, w, rmV);
				assert(ret.size() > 3);
				for (int nx : G.g[w]) {
					if (mk[nx]) rmV[nx] = 0;
				}
				for (int nx : G.g[v]) {
					mk[nx] = 0;
				}
				break;
			}
		}
		if (ret.size() > 0) break;
	}
	if (ret.size() > 0) assert(ret.size() > 3);
	for (int x : order) {
		index[x] = 0;
		f[x] = 0;
	}
	return ret;
}
// returns a chordless cycle in G in component of sv in O(n) time
vector<int> getChordlessCycle(const Graph& G, int sv) {
	assert(G.g[sv].size() > 0);
	vector<int> rmV(G.n);
	return getChordlessCycle(G, sv, rmV);
}

}