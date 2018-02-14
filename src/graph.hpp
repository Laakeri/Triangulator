#ifndef MINFILL_GRAPH_HPP
#define MINFILL_GRAPH_HPP

#include <vector>
#include <ostream>
#include <cstdint>
#include <cassert>
#include <set>

class Graph {
private:
	std::vector<int> vMap;
	std::vector<std::vector<char> > mat;
	bool DENSE;
	void initMatrix();
	bool denseG(int nn, int mm) const;
public:
	std::vector<std::vector<int> > g;
	int n;
	int m;
	Graph(int n_, int m_);
	Graph(int n_, bool dense = false);
	Graph(const std::set<std::pair<int, int> >& es);
	Graph(const std::vector<std::pair<int, int> >& es);
	bool isDense() const;
	void addEdge(int a, int b);
	void popEdge(int a, int b);
	void addEdges(const std::vector<std::pair<int, int> >& es);
	void remEdge(int a, int b);
	void remIncident(int x);
	void remIncident(int x, std::vector<std::pair<int, int> >& remd);
	bool hasEdge(int a, int b) const;
	void print(std::ostream& out) const;
	uint64_t hash() const;
	int mapBack(int v) const;
	std::vector<int> mapBack(std::vector<int> v) const;
	std::pair<int, int> mapBack(std::pair<int, int> e) const;
	std::vector<std::pair<int, int> > mapBack(std::vector<std::pair<int, int> > es) const;
	std::pair<std::vector<std::vector<int> >, std::vector<std::pair<int, int> > > mapBack(std::pair<std::vector<std::vector<int> >, std::vector<std::pair<int, int> > > treeDecomposition) const;
	void dfs1(int x, std::vector<int>& u) const;
	bool isConnected() const;
	bool isConnectedOrIsolated() const;
	bool hasCycle() const;
	int mapInto(int v) const;
	std::vector<std::vector<int> > mapInto(const std::vector<std::vector<int> >& vs) const;
	template <typename T>
	std::vector<T> mapTableInto(const std::vector<T>& table) const {
		std::vector<T> res;
		for (int x : vMap) {
			assert(x >= 0 && x < (int)table.size());
			res.push_back(table[x]);
		}
		return res;
	}
};

#endif