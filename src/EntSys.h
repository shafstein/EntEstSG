#pragma once
#include "EstEntSG.h"
#include "DynSys4EntSys.h"
#include "Maximization.h"
using namespace std;
using namespace arma;

class BasisS {
public:
	int dim;
	BasisS(int _dim);
	cube TheBasis;
	mat E(int nr) { return TheBasis.slice(nr); }
};

umat MulInds(int MaxDegree, int dim);

class EntSys {
public:
	DynSys4EntSys *pSys;
	umat r_a; // the rows of r_a are multiindices
	int dim, NrThreads, PolyDim;
	double Min4Norm;
	mat p, ip, sp, isp;
	cube F, X;
	xGrid xG;
	vec a;
	BasisS Sbas;
	bool RefinedSearch;
	void SetP(const mat &_P);
	void Set_a(const vec &_a) { a = _a; };
	void UpdateF_X(void);
	double TrNormSQ(const mat &s2);
	void Norms1s2(vec &s1, mat &s2);
	virtual tuple<int, double> k0EntEst(const vec &lambda, double cfac) = 0;
	// returns < s1n,s2n,k0,EntEst > with s1,s2 dummies is OnlyEntEst
	virtual tuple<vec, mat, int, double> riem_subg(const vec &cx, bool OnlyEntEst = false) = 0;
	virtual mat Bmat(const mat &A) = 0;
	// pass a pointer to the system and a grid with pSys->CoordTrand(Grid) forward invariant
	EntSys(DynSys4EntSys *_pSys, umat _r_a, const uvec &resol, bool _RefinedSearch = true, int NrThreads = -1, double _Min4Norm = 1e-5);
	// returns the maximum and pSys
	tuple<double, vec> FindMaximum(void);
	virtual double abl(const vec &ox) = 0;
	void StepForward(double t, const vec &s1, const mat &s2);
	virtual vec s1_vec(const vec &cx) = 0;
};

class EntSysCONT : public EntSys	// continuous time systems
{
public:
	tuple<int, double> k0EntEst(const vec &lambda, double cfac);
	mat Bmat(const mat &A);
	tuple<vec, mat, int, double> riem_subg(const vec &cx, bool OnlyEntEst = false);
	double abl(const vec &ox);
	vec s1_vec(const vec &cx);
	EntSysCONT(DynSys4EntSys *_pSys, umat _r_a, const uvec &resol, bool _RefinedSearch = true, int _NrThreads = -1)
		: EntSys(_pSys, _r_a, resol, _RefinedSearch, _NrThreads) {};
};

class EntSysDISC : public EntSys	// discrete time systems
{
public:
	tuple<int, double> k0EntEst(const vec &lambda, double cfac);
	mat Bmat(const mat &A);
	tuple<vec, mat, int, double> riem_subg(const vec &cx, bool OnlyEntEst = false);
	double abl(const vec &ox);
	EntSysDISC(DynSys4EntSys *_pSys, umat _r_a, const uvec &resol, bool _RefinedSearch = true, int _NrThreads = -1)
		: EntSys(_pSys, _r_a, resol, _RefinedSearch, _NrThreads) {};
	vec s1_vec(const vec &cx);
};
