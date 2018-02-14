#ifndef MINFILL_IO_HPP
#define MINFILL_IO_HPP

#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>

class IO {
private:
	std::map<std::string, int> vertices;
	std::vector<std::string> verticeLabels;
	std::string vertexId(int id) const;
	int paceVertices;
	bool paceTreeWidth;
	int nTokens(std::string s) const;
public:
	int readInput(std::istream& in, std::set<std::pair<int, int> >& edges);
	void printFill(std::ostream& out, const std::vector<std::pair<int, int> >& edges) const;
	void printTreeDecomposition(std::ostream& out, const std::pair<std::vector<std::vector<int> >, std::vector<std::pair<int, int> > > treeDecomposition);
};

#endif