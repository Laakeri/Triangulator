#ifndef MINFILL_GLOBAL_HPP
#define MINFILL_GLOBAL_HPP

#include <vector>
#include <cstdint>
#include "timer.hpp"

const int maxN = 50000;
const int maxN2 = 1000000;
const std::vector<std::pair<int, int> > FAILSOLUTION = {{-1, -1}};

const int64_t maxBWidth = 1e18;

const int PROBLEM_MINFILL = 1;
const int PROBLEM_TREEWIDTH = 2;
const int PROBLEM_BAYES = 3;
const int PROBLEM_HYPER = 4;
const int PROBLEM_DECIDE = 5;
const int PROBLEM_FRACHYPER = 6;

#endif