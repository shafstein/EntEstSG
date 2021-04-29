#include "Maximization.h"

using namespace std;
using namespace arma;

ZGrid::ZGrid(const ivec &_L, const ivec &_U, const int _Gdim) : L(_L), U(_U), Gdim(_Gdim) {
	assert(L.n_rows == U.n_rows && U.n_rows == Gdim);
	bint OldEndIndex;
	if (min(U - L) < 0) { // empty grid
		EndIndex = 0;
	} else { // grid is not empty, calculate number of elements
		EndIndex = 1;
		for (int i = 0; i < Gdim; i++) {
			OldEndIndex = EndIndex;
			EndIndex *= U(i) - L(i) + 1;
			assert(OldEndIndex <= EndIndex); // check for overflow
		}
	}
}

ZGrid::~ZGrid() {}

bint ZGrid::V2I(ivec v) const {
	assert(InGrid(v));
	bint i, Mul, Index;
	v -= L;
	for (i = 0, Mul = 1, Index = 0; i < Gdim; i++) {
		Index += v(i) * Mul;
		Mul *= U(i) - L(i) + 1;
	}
	return Index;
}

vector<bint> ZGrid::V2I(vector<ivec> v) const {
	vector<bint> iv;
	for (auto val = v.begin(); val < v.end(); val++) {
		iv.push_back(V2I(*val));
	}
	return iv;
}

ivec ZGrid::I2V(bint Index) const {
	if (!(0 <= Index && Index < EndIndex)) {
		cout << Index << " y " << EndIndex << endl;
		cout << (0 <= Index) << " " << (Index < EndIndex) << endl;
	}
	assert(0 <= Index && Index < EndIndex);
	ivec v(Gdim);
	for (int i = 0; i < Gdim; i++) {
		v(i) = Index % (U(i) - L(i) + 1);
		Index /= U(i) - L(i) + 1;
	}
	return v += L;
}

vec ZGrid::I2vec(bint Index) const {
	return conv_to<vec>::from(I2V(Index));
}

bool ZGrid::InGrid(const ivec &v) const {
	assert(v.n_rows == Gdim);
	return min(v - L) >= 0 && max(v - U) <= 0;
}

// xGrid

xGrid::xGrid(const vec &_xL, const vec &_xU, const uvec &iU, const int _Gdim) : xL(_xL), xU(_xU), G(zeros<ivec>(_Gdim), conv_to<ivec>::from(iU), _Gdim), Gdim(_Gdim) {
	assert(xL.n_rows == xU.n_rows && xU.n_rows == Gdim);
	assert(min(xU - xL) >= 0);
	assert(min(iU) >= 1);
	hv = (xU - xL) / conv_to<vec>::from(iU);
}

bool xGrid::InArea(const arma::vec &x) const {
	return min(x - xL) >= 0.0 && min(xU - x) >= 0.0;
}

vec xGrid::I2vec(bint Index) const {
	return xL + conv_to<vec>::from(G.I2V(Index)) % hv;
}

double MaxInXGridParallel(const xGrid &xG, bint &PosMax, function<double(const arma::vec &)> ObjFunc, int NrThreads) { // ox=I2V(PosMax), cx=CoordTrans(ox)
	assert(0 < NrThreads);
	bint NrProb = xG.NrPoints();
	vector<bint> split(NrThreads + 1);
	vector<bint> vPosMax(NrThreads);
	vec vMaxVal(NrThreads);
	for (int j = 0; j < NrThreads; j++) {
		split[j] = j * NrProb / NrThreads;
	}
	split[NrThreads] = NrProb;
	function<void(int)> SearchInterval = [&vMaxVal, &split, &xG, &ObjFunc, &vPosMax](int j) {
		double TvMaxVal = -Infinity;
		bint TvPosMax = -1;
		for (bint i = split[j]; i < split[j + 1]; i++) {
			vec x = xG.I2vec(i);
			double val = ObjFunc(x);
			if (val > TvMaxVal) {
				TvMaxVal = val;
				TvPosMax = i;
			}
		}
		vMaxVal[j] = TvMaxVal;
		vPosMax[j] = TvPosMax;
	};

	vector<thread> threads(NrThreads);
	for (int j = 0; j < NrThreads; j++) {
		threads[j] = thread(SearchInterval, j);
	}
	for (int j = 0; j < NrThreads; j++) {
		threads[j].join();
	}
	int MaxInd = vMaxVal.index_max();
	PosMax = vPosMax[MaxInd];
	return vMaxVal[MaxInd];
}