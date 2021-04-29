#include "EntSys.h"

double multind(const vec &x, const uvec &ind) {
	double ret = 1.0;
	for (int i = 0; i < ind.n_elem; i++) {
		if (ind(i) > 0) {
			ret *= pow(x(i), ind(i));
		}
	}
	return ret;
}

BasisS::BasisS(int _dim) : dim(_dim) {
	TheBasis = cube(dim, dim, (dim * (dim + 1)) / 2);
	int m = 0;
	for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			TheBasis.slice(m) = zeros<mat>(dim, dim);
			TheBasis(i, j, m) = TheBasis(j, i, m) = 1;
			m++;
		}
	}
}

umat MulInds(int MaxDegree, int dim) { // make all multiindices >0 up to MaxDegree
	vector<umat> Indices;
	function<void(uvec, int, int)> f = [&](uvec uv, int i, int Nr) -> void {
		if (Nr == MaxDegree || i == dim)
			return;
		uv(i) += 1;
		Indices.push_back(uv);
		f(uv, i, Nr + 1);
		uv(i) -= 1;
		f(uv, i + 1, Nr);
	};
	f(zeros<uvec>(dim), 0, 0);
	int pos = 0;
	int NrIndices = Indices.size();
	umat ret(NrIndices, dim);
	for (int deg = 1; deg <= MaxDegree; deg++) {
		for (int j = 0; j < NrIndices; j++) {
			if (accu(Indices[j]) == deg) {
				ret.row(pos++) = Indices[j].t();
			}
		}
	}
	return ret;
}

// EntSys
EntSys::EntSys(DynSys4EntSys *_pSys, umat _r_a, const uvec &resol, bool _RefinedSearch, int _NrThreads, double _Min4Norm)
	: pSys(_pSys), dim(_pSys->dim), r_a(_r_a), xG(zeros<vec>(dim), ones<vec>(dim), resol, dim), Sbas(_pSys->dim), RefinedSearch(_RefinedSearch), NrThreads(_NrThreads), Min4Norm(_Min4Norm) {
	SetP(eye<mat>(dim, dim));
	PolyDim = r_a.n_rows;
	a = zeros<vec>(PolyDim);
	F = cube(dim, dim, (dim * (dim + 1)) / 2);
	X = cube(dim, dim, (dim * (dim + 1)) / 2);
	if (NrThreads == -1) {
		NrThreads = thread::hardware_concurrency();
	}
	NrThreads = max(min(bint(NrThreads), xG.NrPoints() / 100), bint(1)); // if NrThreads=0
}

tuple<double, vec> EntSys::FindMaximum(void) {
	double MaxVal;
	bint PosMax;
	MaxVal = MaxInXGridParallel(xG, PosMax, [&](const vec &ox) -> double { return abl(ox); }, NrThreads);
	if (RefinedSearch) {
		vec sx = xG.I2vec(PosMax);
		// make sure not to search outside of the area
		vec epssmall = 0.5 * xG.hv;
		mat lbm(dim, 2);
		lbm.col(0) = sx - epssmall;
		lbm.col(1) = xG.xL;
		mat ubm(dim, 2);
		ubm.col(0) = sx + epssmall;
		ubm.col(1) = xG.xU;
		uvec resolxG2 = conv_to<uvec>::from(xG.G.U);
		xGrid xG2(max(lbm, 1), min(ubm, 1), resolxG2, dim);
		MaxVal = MaxInXGridParallel(xG2, PosMax, [&](const vec &ox) -> double { return abl(ox); }, NrThreads);
		return tuple<double, vec>{MaxVal, pSys->CoordTrans(xG2.I2vec(PosMax))};
	}
	else {
		return tuple<double, vec>{MaxVal, pSys->CoordTrans(xG.I2vec(PosMax))};
	}
}

void EntSys::SetP(const mat &_p) {
	p = _p;
	ip = inv(p);
	sp = sqrtmat_sympd(p);
	isp = inv(sp);
}

void EntSys::UpdateF_X(void) {
	F = X = zeros<cube>(dim, dim, (dim * (dim + 1)) / 2);
	F.slice(0) = Sbas.E(0) / sqrt(trace(ip * Sbas.E(0) * ip * Sbas.E(0)));
	for (int i = 1; i < (dim * (dim + 1)) / 2; i++) {
		F.slice(i) = Sbas.E(i);
		for (int j = 0; j < i; j++) {
			F.slice(i) -= trace(ip * F.slice(j) * ip * Sbas.E(i)) * F.slice(j);
		}
		F.slice(i) = F.slice(i) / sqrt(trace(ip * F.slice(i) * ip * F.slice(i)));
	}
	for (int i = 0; i < (dim * (dim + 1)) / 2; i++) {
		X.slice(i) = syl(sp, sp.t(), -F.slice(i));
	}
}

void EntSys::Norms1s2(vec &s1, mat &s2) {
	double s1norm = norm(s1);
	double s2norm = sqrt(abs(trace(ip * s2 * ip * s2))); // abs for safety
	if (s1norm > Min4Norm) {
		s1 = s1 / s1norm;
	}
	if (s2norm > Min4Norm) {
		s2 = s2 / s2norm;
	}
}

void EntSys::StepForward(double t, const vec &s1, const mat &s2) {
	a -= t * s1;
	mat NewP = sp * expmat(t * isp * (-s2) * isp) * sp;
	SetP(NewP);
}

// EntSysCONT
mat EntSysCONT::Bmat(const mat &A) {
	return sp * A * isp + isp * A.t() * sp;
}

tuple<int, double> EntSysCONT::k0EntEst(const vec &lambda, double cfac) {
	int k0 = 0;
	double EntEst = 0.0;
	for (int i = 0; i < dim; i++) {
		double val = lambda(i) + cfac;
		if (val > 0.0) {
			k0++;
			EntEst += val;
		}
	}
	EntEst /= 2.0 * log(2);
	return tuple<int, double>{k0, EntEst};
}

double EntSysCONT::abl(const vec &ox) {
	vec cx = pSys->CoordTrans(ox);
	double cfac = dot(a, s1_vec(cx));
	vec lambda = eig_sym(Bmat(pSys->Amat(cx)));
	double retval = 0.0;
	for (int i = 0; i < dim; i++) {
		double val = cfac + lambda(i);
		if (val > 0.0) {
			retval += val;
		}
	}
	return retval;
}

vec EntSysCONT::s1_vec(const vec &cx) {
	vec fcx = pSys->f(cx);
	vec ret = zeros<vec>(PolyDim);
	for (int i = 0; i < PolyDim; i++) {
		uvec m_row_i = r_a.row(i).t();
		for (int j = 0; j < dim; j++) {
			if (m_row_i(j) > 0) {
				uvec m = m_row_i;
				m(j) -= 1;
				ret(i) += m_row_i(j) * multind(cx, m) * fcx(j);
			}
		}
	}
	return ret;
}

tuple<vec, mat, int, double> EntSysCONT::riem_subg(const vec &cx, bool OnlyEntEst) {
	mat V;
	vec lambda;
	mat A = pSys->Amat(cx);
	eig_sym(lambda, V, Bmat(A));
	// change from ascending to decending order in lambda and adapt V
	mat rev = zeros<mat>(dim, dim);
	for (int i = 0; i < dim; i++) {
		rev(dim - 1 - i, i) = 1.0;
	}
	lambda = rev * lambda;
	V = V * rev;

	vec s1 = s1_vec(cx);
	double cfac = dot(a, s1);
	int k0;
	double EntEst;
	tie(k0, EntEst) = k0EntEst(lambda, cfac);
	if (OnlyEntEst || k0 == 0) {
		return tuple<vec, mat, int, double>{zeros<vec>(PolyDim), zeros<mat>(dim, dim), k0, EntEst};
	}
	s1 *= k0;

	// compute s2
	mat D = zeros<mat>(dim, dim);
	for (int i = 0; i < k0; i++) {
		D(i, i) = 1.0;
	}
	mat S = V * D * V.t();
	mat s2 = zeros<mat>(dim, dim);
	UpdateF_X();
	for (int i = 0; i < (dim * (dim + 1)) / 2; i++) {
		s2 += trace(S.t() * (X.slice(i) * A * isp - sp * A * isp * X.slice(i) * isp - isp * X.slice(i) * isp * A.t() * sp + isp * A.t() * X.slice(i))) * F.slice(i);
	}
	return tuple<vec, mat, int, double>{s1, s2, k0, EntEst};
}

// EntSysDISC
mat EntSysDISC::Bmat(const mat &A) {
	return sp * A * isp;
}
tuple<int, double> EntSysDISC::k0EntEst(const vec &alpha, double cfac) {
	int k0 = 0;
	double EntEst = 0.0;
	for (int i = 0; i < dim; i++) {
		double val = log(alpha(i)) + cfac;
		if (val > 0.0) {
			k0++;
			EntEst += val;
		}
	}
	EntEst /= log(2);
	return tuple<int, double>{k0, EntEst};
}
double EntSysDISC::abl(const vec &ox) {
	vec cx = pSys->CoordTrans(ox);
	double cfac = 0.5 * dot(a, s1_vec(cx));
	vec alpha = svd(Bmat(pSys->Amat(cx)));
	double retval = 0.0;
	for (int i = 0; i < dim; i++) {
		double val = cfac + log(alpha(i));
		if (val > 0.0) {
			retval += val;
		}
	}
	return retval;
}

vec EntSysDISC::s1_vec(const vec &cx) {
	vec fcx = pSys->f(cx);
	vec ret(PolyDim);
	for (int i = 0; i < PolyDim; i++) {
		ret(i) = multind(fcx, r_a.row(i).t()) - multind(cx, r_a.row(i).t());
	}
	return ret;
}

tuple<vec, mat, int, double> EntSysDISC::riem_subg(const vec &cx, bool OnlyEntEst) {
	mat U, V;
	vec alpha;
	mat A = pSys->Amat(cx);
	svd(U, alpha, V, Bmat(A));
	vec s1 = s1_vec(cx);
	double cfac = 0.5 * dot(a, s1);
	int k0;
	double EntEst;
	tie(k0, EntEst) = k0EntEst(alpha, cfac);
	if (OnlyEntEst || k0 == 0) {
		return tuple<vec, mat, int, double>{zeros<vec>(PolyDim), zeros<mat>(dim, dim), k0, EntEst};
	}
	s1 *= 0.5 * k0;

	// compute s2
	UpdateF_X();
	mat D = zeros<mat>(dim, dim);
	for (int i = 0; i < k0; i++) {
		D(i, i) = 1 / alpha(i);
	}
	mat S = U * D * V.t();
	mat s2 = zeros<mat>(dim, dim);
	UpdateF_X();
	for (int i = 0; i < (dim * (dim + 1)) / 2; i++) {
		s2 += trace(S.t() * (X.slice(i) * A * isp - sp * A * isp * X.slice(i) * isp)) * F.slice(i);
	}
	return tuple<vec, mat, int, double>{s1, s2, k0, EntEst};
}