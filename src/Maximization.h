#pragma once
#include"EstEntSG.h"

using namespace std;
using namespace arma;



// for L,U in Z^n all z in Z^n fulfilling L <= z <= U component-wise
class ZGrid {
public:
	const int Gdim;
	arma::ivec L, U;
	bint EndIndex;
	ZGrid(const arma::ivec &_L, const arma::ivec &_U, const int _Gdim);
	~ZGrid();
	bint V2I(arma::ivec) const;
	vector<bint> V2I(vector<arma::ivec>) const;
	arma::ivec I2V(bint) const;
	arma::vec I2vec(bint) const;
	bool InGrid(const arma::ivec &v) const;
	bint NrPoints(void) const { return EndIndex; };
	bool IsEmpty(void) const { return EndIndex == 0; }
};


class xGrid {
public:
	const int Gdim;
	ZGrid G;
	arma::vec xL, xU, hv;
	bool InArea(const arma::vec &x) const;
	arma::vec I2vec(bint Index) const;
	bint NrPoints(void) const { return G.NrPoints(); }
	bint NrCloseGPoint(arma::vec x);
	xGrid(const arma::vec &_xL, const arma::vec &_xU, const arma::uvec &iU, const int _Gdim);
	~xGrid() {};
};

double MaxInXGridParallel(const xGrid &xG, bint &PosMax, function<double(const arma::vec &)> f, int NrThreads);

