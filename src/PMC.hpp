#ifndef MINFILL_PMC_HPP
#define MINFILL_PMC_HPP

#include "graph.hpp"

#include <vector>
#include <cstdint>

namespace PMC {
	struct valueT {
		static const int TYPE_INT = 1;
		static const int TYPE_DOUBLE  = 2;
		int64_t valInt;
		long double valDouble;
		int TYPE;
		valueT();
		valueT(long double x);
		valueT(int64_t x);
		void setV(long double x);
		void setV(int64_t x);
		bool operator<(const valueT& o) const;
		bool operator>(const valueT& o) const;
		bool operator<=(const valueT& o) const;
		bool operator>=(const valueT& o) const;
		valueT operator+(const valueT& o) const;
	};
	std::vector<std::pair<int, int> > solve(const Graph& G, valueT lowerBound, valueT upperBound, int PROBLEM);
	bool hasMoreMinSeps(const Graph& G, int64_t limit);
	void initTableSz(std::vector<int64_t> ts);
	void initVertMap(std::vector<int> vertMap);
	std::pair<int, int> getMinSepPMC();
	
	void setPrintProgress(bool val);
};

#endif