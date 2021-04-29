#pragma once
#include "EstEntSG.h"

using namespace std;
using namespace arma;

class DynSys4EntSys {
public:
	int dim;
	virtual vec f(const vec &cx) = 0;
	virtual vec CoordTrans(const vec &ox) =0;
	virtual mat Amat(const vec &cx) = 0;	  // Df(cx)
	DynSys4EntSys(int _dim) : dim(_dim) {};
};


// concrete systems

class Lorenz : public DynSys4EntSys {
public:
	double sigma, rho, beta, K_Rad;
	Lorenz(double _sigma = 10.0, double _rho = 28.0, double _beta = 8.0 / 3.0);
	vec f(const vec &cx);
	vec CoordTrans(const vec &ox);
	mat Amat(const vec &cx);
};

class BouncingBall : public DynSys4EntSys {
public:
	double delta, gamma;
	BouncingBall(double _delta = 2.0, double _gamma = 0.1)
		: delta(_delta), gamma(_gamma), DynSys4EntSys(2) {};
	vec f(const vec &cx);
	vec CoordTrans(const vec &ox); 
	mat Amat(const vec &cx);
};

class Henon : public DynSys4EntSys {
public:
	double TransLambda;
	vec Transq00, Transl;
	mat TranssMat;
	Henon(void);
	vec f(const vec &cx);
	vec CoordTrans(const vec &ox);
	mat Amat(const vec &cx);
};
