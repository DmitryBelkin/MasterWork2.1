#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <set>
#include <fstream>

using namespace std;

const string outputPrefix = "../resources/output/";
const string solverInfoFilename = outputPrefix + "solverInfo.txt";

class MCG
{
public:
	unsigned int         n, m_maxiter;
	double               m_eps;
	vector < int >       ia, ja;
	vector < double >    xtch, r, z, temp, L, U, LUdi, t, f, ggl, ggu, di;
	vector < set <int> > ig;

	MCG()
		: n(0)
		, m_maxiter(100000)
		, m_eps(1e-09)
	{ 
		ofstream solverInfo(solverInfoFilename, ios_base::trunc);
	};
	~MCG() { };

	void   MCG_LU();
	void   CreateLU();

	void   SetEps(const double eps);
	void   SetMaxiter(const int maxiter);

	double NormVector(const vector <double> &x) const;
	double ScalarProduct(const vector <double> &x, const vector <double> &y) const;

	void   Lx (const vector <double> &ELf, const vector <double> &Df, vector <double> &x, const vector <double> &f) const;
	void   Ux (const vector <double> &EUf,                            vector <double> &x, const vector <double> &f) const;
	void   LTx(const vector <double> &ELf, const vector <double> &Df, vector <double> &x, const vector <double> &f) const;
	void   UTx(const vector <double> &EUf,                            vector <double> &x, const vector <double> &f) const;
	void   MultMatrixOnVector(const vector <double> &EU, const vector <double> &EL, const vector <double> &D, const vector <double> &vect, vector <double> &res) const;

	void WriteSolverInfoInFile   (const int iterations, const double residual) const;
	void WriteSolverInfoInConsole(const int iterations, const double residual) const;
};
