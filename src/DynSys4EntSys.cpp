#include "DynSys4EntSys.h"

// Lorenz
Lorenz::Lorenz(double _sigma, double _rho, double _beta) : sigma(_sigma), rho(_rho), beta(_beta), DynSys4EntSys(3) {
	K_Rad = sqrt(beta / 2.0) * (sigma + rho);
}

vec Lorenz::f(const vec &cx) {
	vec fx(dim);
	fx(0) = sigma * (-cx(0) + cx(1));
	fx(1) = rho * cx(0) - cx(1) - cx(0) * cx(2);
	fx(2) = -beta * cx(2) + cx(0) * cx(1);
	return fx;
}

mat Lorenz::Amat(const vec &cx) {
	return { {-sigma, sigma, 0}, {rho - cx(2), -1, -cx(0)}, {cx(1), cx(0), -beta} };
}

vec Lorenz::CoordTrans(const vec &ox) {
	vec sox = ox % vec{ K_Rad, datum::pi, 2.0 * datum::pi };
	vec cx(dim);
	cx(0) = sox(0) * sin(sox(1)) * cos(sox(2));
	cx(1) = sox(0) * sin(sox(1)) * sin(sox(2));
	cx(2) = sox(0) * cos(sox(1)) + (sigma + rho);
	return cx;
}

// BouncingBall
vec BouncingBall::f(const vec &cx) {
	vec fx(dim);
	fx(0) = cx(0) + cx(1);
	fx(1) = gamma * cx(1) - delta * cos(cx(0) + cx(1));
	return fx;
}

mat BouncingBall::Amat(const vec &cx) {
	return { {1.0, 1.0}, {delta * sin(cx(0) + cx(1)), gamma + delta * sin(cx(0) + cx(1))} };
}

vec BouncingBall::CoordTrans(const vec &ox) {
	double x_min = -delta / (1 - gamma), x_max = 2.0 * datum::pi, y_max = delta / (1.0 - gamma);
	vec cx(dim);
	cx(0) = x_min + ox(0) * (x_max - x_min);
	cx(1) = ox(1) * y_max;
	return cx;
}

// Henon a=1.4, b=0.3
Henon::Henon(void) : DynSys4EntSys(2) {
	vec q00 = { -1.862, 1.96 };
	vec q10 = { -1.484, -2.3333 };
	vec q11 = { 1.743, -0.6533 };
	vec q01 = { 1.848, 0.6267 };
	mat Mq = { {q10(0) - q00(0), q01(0) - q00(0)}, {q10(1) - q00(1), q01(1) - q00(1)} };
	vec cd = solve(Mq, q11 - q00);
	TransLambda = cd(0) + cd(1) - 1;
	Transl = vec{ 1 - cd(1), 1 - cd(0) };
	TranssMat = Mq * diagmat(cd) + q00 * Transl.t();
	Transq00 = q00;
	// (TranssMat * x + Translambda * Transq00)/ (dot(Transl, x)  + Translambda)
}

vec Henon::f(const vec &cx) {
	vec fx(dim);
	fx(0) = 1.4 - pow(cx(0), 2) + 0.3 * cx(1);
	fx(1) = cx(0);
	return fx;
}

mat Henon::Amat(const vec &cx) {
	return { {-2.0 * cx(0), 0.3}, {1.0, 0.0} };
}

vec Henon::CoordTrans(const vec &ox) {
	return (TranssMat * ox + TransLambda * Transq00) / (dot(Transl, ox) + TransLambda);
}