#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

class MCG
{
public:
	int                    n, m_maxiter;
	double                 m_eps;
	vector < int >         ia, ja;
	vector < double >      xtch, r, z, temp, L, U, LUdi, t, f, ggl, ggu, di;
	vector <vector <int> > ig;

	MCG() { };
	~MCG() { };

	void setMaxiter(int maxiter);
	void setEps(double eps);

	void Lx(vector <double> &ELf, vector <double> &Df, vector <double> &x, vector <double> &f);
	void Ux(vector <double> &EUf, vector <double> &x, vector <double> &f);
	void LTx(vector <double> &ELf, vector <double> &Df, vector <double> &x, vector <double> &f);
	void UTx(vector <double> &EUf, vector <double> &x, vector <double> &f);
	double scalarProduct(vector <double> &x, vector <double> &y);
	double normVector(vector <double> &x);
	void multMatrixOnVector(vector <double> &EU, vector <double> &EL, vector <double> &D, vector <double> &vect, vector <double> &res);
	void createLU();
	void MCG_LU();
};

