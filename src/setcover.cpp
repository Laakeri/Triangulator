#include "setcover.hpp"
#include "timer.hpp"
#include "global.hpp"

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>

#ifdef USE_GLPK
#include <glpk.h>
#else
#endif

#define F first
#define S second

using namespace std;

typedef long double ld;

namespace SetCover {

vector<int> cv;
vector<int> opt;

int gm;

int git = 1;
int sel[maxN];
int cnt[maxN];
int plc[maxN];
int is[maxN];
int heur[maxN];
int um[maxN];
pair<int, int> cc[maxN];
vector<int> where[maxN];
vector<int> Gsets[maxN];
vector<int> sets[maxN];

int no[maxN];
int tno[maxN];

int tnoS = 0;

bool hasInit = false;

void init(vector<vector<int> > ss) {
	assert(!hasInit);
	hasInit = true;
	gm = ss.size();
	for (int i = 0; i < gm; i++) {
		Gsets[i] = ss[i];
	}
}

int okCover;
int upperBound;

void go(int i, int n) {
	if (i == n) {
		assert((int)cv.size() <= upperBound);
		upperBound = (int)cv.size();
		if (opt.size() == 0) {
			opt = cv;
		}
		else {
			assert((int)cv.size() < (int)opt.size());
			opt = cv;
		}
		assert(opt.size() > 0);
		return;
	}
	bool f = false;
	for (int x : where[i]) {
		if (sel[x]) {
			f = true;
			break;
		}
	}
	if (!f) {
		int ss = tnoS;
		for (int x : where[i]) {
			if ((int)cv.size() + 1 > upperBound) break;
			if ((int)cv.size() + 1 == (int)opt.size()) break;
			if (opt.size() != 0 && (int)opt.size() <= okCover) break;
			if (no[x]) continue;
			sel[x] = 1;
			cv.push_back(x);
			go(i+1, n);
			sel[x] = 0;
			cv.pop_back();
			no[x] = 1;
			tno[tnoS++] = x;
		}
		while (tnoS > ss) {
			no[tno[tnoS - 1]] = 0;
			tnoS--;
		}
	}
	else {
		go(i+1, n);
	}
}

Timer SCTimer1;
Timer SCTimer2;

double getTime1() {
	return SCTimer1.getTime().count();
}
double getTime2() {
	return SCTimer2.getTime().count();
}

map<vector<int>, vector<int> > hasSol;
map<vector<int>, int> noSol;

vector<int> solve(vector<int> univ, int okCover_, int upperBound_) {
	git++;
	assert(hasInit);
	if (univ.size() == 0) return {};
	SCTimer1.start();
	okCover = okCover_;
	upperBound = upperBound_;
	assert(okCover <= upperBound);
	sort(univ.begin(), univ.end());
	univ.erase(unique(univ.begin(), univ.end()), univ.end());
	
	if (hasSol.count(univ) && (int)hasSol[univ].size() <= okCover) {
		return hasSol[univ];
	}
	if (noSol.count(univ) && noSol[univ] >= upperBound) {
		return {-1};
	}
	
	int n = univ.size();
	
	for (int i = 0; i < n; i++) {
		assert(univ[i] >= 0 && univ[i] < maxN);
		plc[univ[i]] = i;
		is[univ[i]] = git;
	}
	int m = 0;
	for (int i = 0; i < gm; i++) {
		sets[m].clear();
		for (int j = 0; j < (int)Gsets[i].size(); j++) {
			if (is[Gsets[i][j]] == git) {
				sets[m].push_back(plc[Gsets[i][j]]);
			}
		}
		if (sets[m].size() > 0) {
			for (int j = 0; j < (int)sets[m].size(); j++) {
				if (j > 0) assert(sets[m][j] > sets[m][j-1]);
				assert(sets[m][j] >= 0 && sets[m][j] < n);
			}
			m++;
		}
	}
	sort(sets, sets+m);
	auto last = unique(sets, sets+m);
	m = last - sets;
	assert(m <= gm);
	for (int i = 1; i < m; i++) {
		assert(sets[i] > sets[i-1]);
	}
	for (int i = 0; i < n; i++) {
		cnt[i] = 0;
	}
	for (int i = 0; i < m; i++) {
		for (int x : sets[i]) {
			cnt[x]++;
		}
	}
	for (int i = 0; i < n; i++) {
		assert(cnt[i] > 0);
		cc[i] = {cnt[i], i};
	}
	sort(cc, cc + n);
	for (int i = 0; i < n; i++) {
		um[cc[i].S] = i;
		where[i].clear();
	}
	for (int i = 0; i < m; i++) {
		for (int& x : sets[i]) {
			x = um[x];
			where[x].push_back(i);
		}
		sort(sets[i].begin(), sets[i].end());
	}
	for (int i = 0; i < n; i++) {
		for (int x : where[i]) {
			int f = 0;
			for (int j = 0; j < (int)sets[x].size(); j++) {
				if (sets[x][j] == i) {
					f = (int)sets[x].size() - j;
					break;
				}
			}
			assert(f > 0);
			heur[x] = f;
		}
		assert(cc[i].F == (int)where[i].size());
		auto cmp = [&](int a, int b) {
			return heur[a] > heur[b];
		};
		sort(where[i].begin(), where[i].end(), cmp);
	}
	for (int i = 0; i < m; i++) sel[i] = 0;
	opt.clear();
	cv.clear();
	if (hasSol.count(univ) && (int)hasSol[univ].size() <= upperBound) {
		opt = hasSol[univ];
		upperBound = opt.size();
	}
	SCTimer2.start();
	go(0, n);
	SCTimer2.stop();
	if (opt.size() == 0) {
		opt = {-1};
		noSol[univ] = upperBound;
	}
	else {
		if (hasSol.count(univ)) {
			if (opt.size() < hasSol[univ].size()) {
				hasSol[univ] = opt;
			}
		}
		else {
			hasSol[univ] = opt;
		}
	}
	SCTimer1.stop();
	return opt;
}

struct LPSolver {
	const ld eps=1e-9;
	const ld inf=1e30;
	int m, n;
	vector<int> N, B;
	vector<vector<ld> > D;
	LPSolver(vector<vector<ld> >& A, vector<ld>& b, vector<ld>& c) :
	m(b.size()), n(c.size()), N(n + 1), B(m), D(m + 2, vector<ld>(n + 2)) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) D[i][j] = A[i][j];}
		for (int i = 0; i < m; i++) { 
			B[i] = n + i; D[i][n] = -1; D[i][n + 1] = b[i];}
		for (int j = 0; j < n; j++) { N[j] = j; D[m][j] = -c[j]; }
		N[n] = -1; D[m + 1][n] = 1;
	}
	void Pivot(int r, int s) {
		ld inv = 1.0 / D[r][s];
		for (int i = 0; i < m + 2; i++) if (i != r)
			for (int j = 0; j < n + 2; j++) if (j != s)
				D[i][j] -= D[r][j] * D[i][s] * inv;
		for (int j = 0; j < n + 2; j++) if (j != s) D[r][j] *= inv;
		for (int i = 0; i < m + 2; i++) if (i != r) D[i][s] *= -inv;
		D[r][s] = inv;
		swap(B[r], N[s]);
	}
	bool Simplex(int phase) {
		int x = phase == 1 ? m + 1 : m;
		while (true) {
			int s = -1;
			for (int j = 0; j <= n; j++) {
				if (phase == 2 && N[j] == -1) continue;
				if (s==-1||D[x][j]<D[x][s]||(D[x][j]==D[x][s]&&N[j]<N[s])) s=j;
			}
			if (D[x][s] > -eps) return true;
			int r = -1;
			for (int i = 0; i < m; i++) {
				if (D[i][s] < eps) continue;
				if (r == -1 || D[i][n + 1] / D[i][s] < D[r][n + 1] / D[r][s] ||
					((D[i][n + 1]/D[i][s])==(D[r][n+1]/D[r][s])&&B[i]<B[r]))r=i;
			}
			if (r == -1) return false;
			Pivot(r, s);
		}
	}
	ld Solve(vector<ld>& x) {
		int r = 0;
		for (int i = 1; i < m; i++) if (D[i][n + 1] < D[r][n + 1]) r = i;
		if (D[r][n + 1] < -eps) {
			Pivot(r, n);
			if (!Simplex(1) || D[m + 1][n + 1] < -eps) return -inf;
			for (int i = 0; i < m; i++) if (B[i] == -1) {
				int s = -1;
				for (int j = 0; j <= n; j++)
				if (s==-1||D[i][j]<D[i][s]||(D[i][j]==D[i][s]&&N[j]<N[s])) s=j;
				Pivot(i, s);
			}
		}
		if (!Simplex(2)) return inf;
		x = vector<ld>(n);
		for (int i = 0; i < m; i++) if (B[i] < n) x[B[i]] = D[i][n + 1];
		return D[m][n + 1];
	}
};

map<vector<int>, ld> fracSol;

int glpTempId[maxN];
double glpTempVal[maxN];

ld solveFrac(vector<int> univ) {
	git++;
	assert(hasInit);
	if (univ.size() == 0) return {};
	sort(univ.begin(), univ.end());
	univ.erase(unique(univ.begin(), univ.end()), univ.end());
	
	if (fracSol.count(univ)) {
		return fracSol[univ];
	}
	
	int n = univ.size();
	
	for (int i = 0; i < n; i++) {
		assert(univ[i] >= 0 && univ[i] < maxN);
		plc[univ[i]] = i;
		is[univ[i]] = git;
	}
	int m = 0;
	for (int i = 0; i < gm; i++) {
		sets[m].clear();
		for (int j = 0; j < (int)Gsets[i].size(); j++) {
			if (is[Gsets[i][j]] == git) {
				sets[m].push_back(plc[Gsets[i][j]]);
			}
		}
		if (sets[m].size() > 0) {
			for (int j = 0; j < (int)sets[m].size(); j++) {
				if (j > 0) assert(sets[m][j] > sets[m][j-1]);
				assert(sets[m][j] >= 0 && sets[m][j] < n);
			}
			m++;
		}
	}
	sort(sets, sets+m);
	auto last = unique(sets, sets+m);
	m = last - sets;
	assert(m <= gm);
	for (int i = 1; i < m; i++) {
		assert(sets[i] > sets[i-1]);
	}
	for (int i = 0; i < n; i++) {
		where[i].clear();
	}
	for (int i = 0; i < m; i++) {
		for (int x : sets[i]) {
			where[x].push_back(i);
		}
	}
#ifdef USE_GLPK
	glp_prob* lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MIN);
	assert(glp_add_rows(lp, n) == 1);
	assert(glp_add_cols(lp, m) == 1);
	for (int i = 0; i < n; i++) {
		glp_set_row_bnds(lp, i+1, GLP_LO, 1, 0);
	}
	for (int i = 0; i < m; i++) {
		glp_set_col_bnds(lp, i+1, GLP_DB, 0, 1);
		glp_set_obj_coef(lp, i+1, 1);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < (int)where[i].size(); j++) {
			glpTempId[j+1] = where[i][j]+1;
			glpTempVal[j+1] = 1;
		}
		glp_set_mat_row(lp, i+1, (int)where[i].size(), glpTempId, glpTempVal);
	}
	glp_smcp param;
	glp_init_smcp(&param);
	param.msg_lev = GLP_MSG_OFF;
	glp_simplex(lp, &param);
	ld val = glp_get_obj_val(lp);
#else
	vector<vector<ld> > A(n);
	vector<ld> b(n);
	for (int i = 0; i < n; i++) {
		A[i].resize(m);
		for (int x : where[i]) {
			A[i][x] = -1;
		}
		b[i] = -1;
	}
	vector<ld> c(m);
	for (int i = 0; i < m; i++) {
		c[i] = -1;
	}
	LPSolver lps(A, b, c);
	vector<ld> sol;
	ld val = -lps.Solve(sol);
#endif
	fracSol[univ] = val;
	return val;
}
}