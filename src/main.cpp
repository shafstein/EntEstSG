// main.cpp : Defines the entry point for the application.
//

#include "EntSys.h"
using namespace std;
using namespace arma;

int main(void) {
	int NrThreads = 60;
	wall_clock timer;
	timer.tic();
	
		Lorenz sys;
		umat r_a = MulInds(2, sys.dim);
		uvec resol = { 500, 50, 100 };
		EntSysCONT Esys(&sys, r_a, resol, true, NrThreads);
		int N = 4000;			   // no. iterations
		double ta = 2.0, tb = 0.0; // t_j=ta/(j+tb)
		ofstream fout1("./results/LorenzEnt.txt");
		ofstream fout2("./results/LorenzRes.txt");
	
/*
	BouncingBall sys;
	umat r_a = MulInds(3, sys.dim);
	uvec resol = { 1000,1000 };
	EntSysDISC Esys(&sys, r_a, resol, true, NrThreads);
	int N = 200;	// no. iterations
	double ta = 1.0, tb = 0.0;  // t_j=ta/(j+tb)
	ofstream fout1("./results/BBEnt.txt");
	ofstream fout2("./results/BBRes.txt");
*/

/*
	Henon sys;
	umat r_a = MulInds(3, sys.dim);
	uvec resol = { 1000,1000 };
	EntSysDISC Esys(&sys, r_a, resol, true, NrThreads);
	int N = 1000;	// no. iterations
	double ta = 4.0, tb = 0.0;  // t_j=ta/(j+tb)
	ofstream fout1("./results/HenonEnt.txt");
	ofstream fout2("./results/HenonRes.txt");
*/
	int k0, best_k;
	double MaxVal, EstEnt, bestEstEnt = Infinity;
	vec cx, s1, best_a;
	mat s2, best_p;
	double t;
	cout.precision(16);
	fout1.precision(16);
	fout2.precision(16);
	timer.tic();
	for (int j = 0; j < N; j++) {
		tie(MaxVal, cx) = Esys.FindMaximum();
		tie(s1, s2, k0, EstEnt) = Esys.riem_subg(cx);
		cout << "itr " << j << ": MaxVal " << MaxVal << ": Ent est " << EstEnt << endl;
		fout1 << EstEnt << endl;
		if (EstEnt < bestEstEnt) {
			bestEstEnt = EstEnt;
			best_p = Esys.p;
			best_a = Esys.a;
			best_k = j;
		}
		t = ta / (j + 1 + tb);
		Esys.Norms1s2(s1, s2);
		Esys.StepForward(t, s1, s2); // step j+1
	}
	cout << N << " iterations computed in " << timer.toc() << " s" << endl;
	fout2 << N << " iterations computed in " << timer.toc() << " sec. : t_j = " << ta << "/(j+" << tb << ")" << endl;
	fout2 << "Best estimate of Restoration Entropy " << bestEstEnt << endl
		<< "obtained in iteration " << best_k << " with " << endl
		<< "a : (pow of x_1, pow of x_2, etc. )  : coefficient " << endl;
	for (int i = 0; i < Esys.PolyDim; i++) {
		fout2 << "( ";
		for (int j = 0; j < Esys.dim; j++) {
			fout2 << Esys.r_a(i, j) << " ";
		}
		fout2 << ") : " << best_a(i) << endl;
	}
	fout2 << endl
		<< "p :" << endl;
	best_p.raw_print(fout2);
}