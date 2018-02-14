#include "graph.hpp"

#include <vector>
#include <ostream>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <iostream>

#define F first
#define S second

using namespace std;

void Graph::initMatrix() {
	mat.resize(n);
	for (int i = 0; i < n; i++) {
		mat[i].resize(n);
	}
}
bool Graph::denseG(int nn, int mm) const {
	return ((int64_t)nn*nn*nn <= (int64_t)mm*mm*20);
}
Graph::Graph(int n_, int m_) {
	n = n_;
	m = 0;
	g.resize(n);
	DENSE = denseG(n, m_);
	if (DENSE) initMatrix();
}
Graph::Graph(int n_, bool dense) {
	n = n_;
	m = 0;
	g.resize(n);
	DENSE = dense;
	if (DENSE) initMatrix();
}
Graph::Graph(const set<pair<int, int> >& es) {
	for (auto e : es) {
		vMap.push_back(e.F);
		vMap.push_back(e.S);
	}
	sort(vMap.begin(), vMap.end());
	vMap.erase(unique(vMap.begin(), vMap.end()), vMap.end());
	n = vMap.size();
	m = 0;
	g.resize(n);
	DENSE = denseG(n, es.size());
	if (DENSE) initMatrix();
	for (auto e : es) {
		e.F = lower_bound(vMap.begin(), vMap.end(), e.F) - vMap.begin();
		e.S = lower_bound(vMap.begin(), vMap.end(), e.S) - vMap.begin();
		addEdge(e.F, e.S);
	}
}
Graph::Graph(const vector<pair<int, int> >& es) {
	for (auto e : es) {
		vMap.push_back(e.F);
		vMap.push_back(e.S);
	}
	sort(vMap.begin(), vMap.end());
	vMap.erase(unique(vMap.begin(), vMap.end()), vMap.end());
	n = vMap.size();
	m = 0;
	g.resize(n);
	DENSE = denseG(n, es.size());
	if (DENSE) initMatrix();
	for (auto e : es) {
		e.F = lower_bound(vMap.begin(), vMap.end(), e.F) - vMap.begin();
		e.S = lower_bound(vMap.begin(), vMap.end(), e.S) - vMap.begin();
		addEdge(e.F, e.S);
	}
}
bool Graph::isDense() const {
	return DENSE;
}
void Graph::addEdge(int a, int b) {
	if (DENSE) {
		mat[a][b] = 1;
		mat[b][a] = 1;
	}
	g[a].push_back(b);
	g[b].push_back(a);
	m++;
}
void Graph::popEdge(int a, int b) {
	if (DENSE) {
		mat[a][b] = 0;
		mat[b][a] = 0;
	}
	g[a].pop_back();
	g[b].pop_back();
	m--;
}
void Graph::addEdges(const vector<pair<int, int> >& es) {
	for (auto e : es) {
		addEdge(e.F, e.S);
	}
}
void Graph::remEdge(int a, int b) {
	if (DENSE) {
		mat[a][b] = 0;
		mat[b][a] = 0;
	}
	for (size_t i = 0; i < g[a].size(); i++) {
		if (g[a][i] == b) {
			swap(g[a][i], g[a].back());
			g[a].pop_back();
			break;
		}
	}
	for (size_t i = 0; i < g[b].size(); i++) {
		if (g[b][i] == a) {
			swap(g[b][i], g[b].back());
			g[b].pop_back();
			break;
		}
	}
	m--;
}
void Graph::remIncident(int x) {
	for (int nx : g[x]) {
		for (size_t i = 0; i < g[nx].size(); i++) {
			if (g[nx][i] == x) {
				swap(g[nx][i], g[nx].back());
				g[nx].pop_back();
				break;
			}
		}
		if (DENSE) {
			mat[x][nx] = 0;
			mat[nx][x] = 0;
		}
		m--;
	}
	g[x].clear();
}
void Graph::remIncident(int x, vector<pair<int, int> >& remd) {
	for (int nx : g[x]) {
		remd.push_back({x, nx});
		for (size_t i = 0; i < g[nx].size(); i++) {
			if (g[nx][i] == x) {
				swap(g[nx][i], g[nx].back());
				g[nx].pop_back();
				break;
			}
		}
		if (DENSE) {
			mat[x][nx] = 0;
			mat[nx][x] = 0;
		}
		m--;
	}
	g[x].clear();
}
bool Graph::hasEdge(int a, int b) const {
	if (DENSE) return mat[a][b] == 1;
	if (g[a].size() < g[b].size()) {
		for (int nx : g[a]) {
			if (nx == b) return true;
		}
		return false;
	}
	else {
		for (int nx : g[b]) {
			if (nx == a) return true;
		}
		return false;
	}
}
void Graph::print(ostream& out) const {
	out<<"n "<<n<<endl;
	int mm = 0;
	for (int i = 0; i < n; i++) {
		mm += (int)g[i].size();
	}
	assert(mm%2 == 0);
	assert(mm/2 == m);
	out<<"m "<<m<<endl;
	for (int i = 0; i < n; i++) {
		for (int nx : g[i]) {
			if (i<nx) out<<i<<" "<<nx<<endl;
		}
	}
}
uint64_t Graph::hash() const {
	uint64_t ha = 0;
	vector<pair<int, int> > es;
	for (int i = 0; i < n; i++) {
		for (int nx : g[i]) {
			if (nx > i) {
				es.push_back({i, nx});
			}
		}
	}
	sort(es.begin(), es.end());
	for (auto e : es) {
		ha *= 1000003ll;
		ha += (uint64_t)e.F;
		ha *= 10007ll;
		ha += (uint64_t)e.S;
	}
	return ha;
}
int Graph::mapBack(int v) const {
	return vMap[v];
}
vector<int> Graph::mapBack(vector<int> v) const {
	for (auto& x : v) {
		x = mapBack(x);
	}
	return v;
}
pair<int, int> Graph::mapBack(pair<int, int> e) const {
	return {mapBack(e.F), mapBack(e.S)};
}
vector<pair<int, int> > Graph::mapBack(vector<pair<int, int> > es) const {
	for (auto& e : es) {
		e = mapBack(e);
	}
	return es;
}
pair<vector<vector<int> >, vector<pair<int, int> > > Graph::mapBack(pair<vector<vector<int> >, vector<pair<int, int> > > treeDecomposition) const {
	for (auto& bag : treeDecomposition.F) {
		bag = mapBack(bag);
	}
	return treeDecomposition;
}
void Graph::dfs1(int x, vector<int>& u) const {
	if (u[x]) return;
	u[x] = 1;
	for (int nx : g[x]) {
		dfs1(nx, u);
	}
}
bool Graph::isConnected() const {
	vector<int> u(n);
	dfs1(0, u);
	for (int i = 0; i < n; i++) {
		if (!u[i]) return false;
	}
	return true;
}
bool Graph::isConnectedOrIsolated() const {
	vector<int> u(n);
	int cc = 0;
	for (int i = 0; i < n; i++) {
		if (u[i] == 0 && g[i].size() > 0) {
			cc++;
			dfs1(i, u);
		}
	}
	return cc <= 1;
}
bool Graph::hasCycle() const {
	vector<int> u(n);
	int cc = 0;
	for (int i = 0; i < n; i++) {
		if (u[i] == 0 && g[i].size() > 0) {
			cc++;
			dfs1(i, u);
		}
	}
	return m > n - cc;
}
int Graph::mapInto(int v) const {
	int id = lower_bound(vMap.begin(), vMap.end(), v) - vMap.begin();
	assert(id < (int)vMap.size());
	assert(vMap[id] == v);
	return id;
}
vector<vector<int> > Graph::mapInto(const vector<vector<int> >& vs) const {
	vector<vector<int> > ret;
	for (auto t : vs) {
		vector<int> r;
		for (int x : t) {
			int id = lower_bound(vMap.begin(), vMap.end(), x) - vMap.begin();
			if (id >= (int)vMap.size()) continue;
			if (vMap[id] != x) continue;
			r.push_back(id);
		}
		sort(r.begin(), r.end());
		if (r.size() > 0) ret.push_back(r);
	}
	sort(ret.begin(), ret.end());
	ret.erase(unique(ret.begin(), ret.end()), ret.end());
	return ret;
}