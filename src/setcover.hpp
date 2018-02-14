#ifndef MINFILL_SETCOVER_HPP
#define MINFILL_SETCOVER_HPP

#include <vector>

namespace SetCover {
	void init(std::vector<std::vector<int> > sets);
	std::vector<int> solve(std::vector<int> univ, int okCover, int upperBound);
	long double solveFrac(std::vector<int> univ);
	double getTime1();
	double getTime2();
}
#endif