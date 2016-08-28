#include "MCG.h"
#include <numeric>

//...........................................................................

void MCG::SetMaxiter(const int maxiter) { m_maxiter = maxiter; }

//...........................................................................

void MCG::SetEps(const double eps) { m_eps = eps; }

//...........................................................................

void MCG::Lx(const vector <double> &ELf, const vector <double> &Df, vector <double> &x, const vector <double> &f) const
{
	vector<double> z(x.size());
	x.swap(z);
	for (unsigned int i = 0; i < n; ++i) x[i] = 0;
	for (unsigned int i = 0; i < n; ++i)
	{
		x[i] = f[i] / Df[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[i] -= ELf[j] * x[ ja[j] ] / Df[i];
	}
}

//...........................................................................

void MCG::Ux(const vector <double> &EUf, vector <double> &x, const vector <double> &f) const
{
	vector<double> z(x.size());
	x.swap(z);
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] += f[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[ ja[j] ] -= EUf[j] * x[i];
	}
}

//...........................................................................

void MCG::LTx(const vector <double> &ELf, const vector <double> &Df, vector <double> &x, const vector <double> &f) const
{
	vector<double> z(x.size());
	x.swap(z);
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = (x[i] + f[i]) / Df[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[ ja[j] ] -= ELf[j] * x[i];
	}
}

//...........................................................................

void MCG::UTx(const vector <double> &EUf, vector <double> &x, const vector <double> &f) const
{
	vector<double> z(x.size());
	x.swap(z);
	for (unsigned int i = 0; i < n; ++i)
	{
		x[i] = f[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[i] -= EUf[j] * x[ ja[j] ];
	}
}

//...........................................................................

double MCG::ScalarProduct(const vector <double> &v1, const vector <double> &v2) const
{
	double temp = 0;
	for (int i = 0; i < n; ++i) temp += v1[i] * v2[i];
	return temp;
	//return inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

//...........................................................................

double MCG::NormVector(const vector <double> &v) const { return sqrt(ScalarProduct(v, v)); }

//...........................................................................

void MCG::MultMatrixOnVector(const vector <double> &EU, const vector <double> &EL, const vector <double> &D, const vector <double> &vect, vector <double> &res) const
{
	vector<double> z(res.size());
	res.swap(z);
	for (unsigned int i = 0; i < n; ++i)
	{
		res[i] = D[i] * vect[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
		{
			res[ ja[j] ] += EU[j] * vect[i];
			res[i]       += EL[j] * vect[ja[j]];
		}
	}
}

//...........................................................................

void MCG::CreateLU()
{
	int size = ia[n], j, kj, ki, j1, k;
	unsigned int i;
	double su, sl;
	L.clear();
	U.clear();
	LUdi.clear();

	L.resize(size);
	U.resize(size);
	LUdi.resize(n);
	for (i = 0; i < n; ++i)
	{
		double sd = 0;
		const int i0 = ia[i];
		const int i1 = ia[i + 1];
		for (k = i0; k < i1; ++k)
		{
			su = 0;
			sl = 0;
			j = ja[k];
			kj = ia[j];
			ki = ia[i];
			j1 = ia[j + 1];
			while ((ki<k) && (kj<j1))
			{
				if (ja[kj] == ja[ki])
				{
					sl += L[ki] * U[kj];
					su += L[kj] * U[ki];
					ki++; kj++;
				}
				else
				{
					if (ja[ki] < ja[kj]) ki++;
					else                 kj++;
				}
			}
			assert(LUdi[j] != 0);
			L[k] =  ggl[k] - sl;
			U[k] = (ggu[k] - su) / LUdi[j];
			sd += L[k] * U[k];
		}
		LUdi[i] = di[i] - sd;
	}

	cout << "5. LU - factorization passed" << endl;
}

//...........................................................................

void MCG::MCG_LU()
{
	unsigned int i, iteration;
	vector <double> r, z, r_, temp1, temp2;
	double t, norm;
	r    .resize(n);
	z    .resize(n);
	r_   .resize(n);
	temp1.resize(n);
	temp2.resize(n);
	CreateLU();
	Lx (L, LUdi, temp1, f);
	LTx(L, LUdi, temp2, temp1);
	MultMatrixOnVector(ggl, ggu, di, temp2, temp1);
	UTx(U, r, temp1);
	z = r_ = r;
	t = NormVector(f);
	cout << "NormVector(r) = " << NormVector(r) << endl;
	cout << "t " << t << endl;
	norm = NormVector(r) / t;
	iteration = 1;
	while (norm > m_eps && iteration < m_maxiter)
	{
		Ux(U, temp1, z);
		MultMatrixOnVector(ggu, ggl, di, temp1, temp2);
		Lx(L, LUdi, temp1, temp2);
		LTx(L, LUdi, temp2, temp1);
		MultMatrixOnVector(ggl, ggu, di, temp2, temp1);
		UTx(U, temp2, temp1);
		const double alpha = ScalarProduct(r, r) / ScalarProduct(temp2, z);
		for (i = 0; i < n; ++i)
		{
			weights[i] += alpha*z[i];
			r[i]       -= alpha*temp2[i];
		}
		const double beta = ScalarProduct(r, r) / ScalarProduct(r_, r_);
		for (i = 0; i < n; ++i)
		{
			z [i]  = r[i] + beta*z[i];
		}
		r_ = r;
		iteration++;
		norm = NormVector(r) / t;

		cout << "\r" << "\t Current residual = " << norm;
	}
	cout << endl;
	Ux(U, temp1, weights);
	weights = temp1;

	cout << "6. Solution obtained" << endl;

	WriteSolverInfoInFile   (iteration, norm);
	WriteSolverInfoInConsole(iteration, norm);
}

//...........................................................................

void MCG::WriteSolverInfoInFile(const int iterations, const double residual) const
{
	ofstream solverInfo(solverInfoFilename, ios::out);
	solverInfo.setf(ios::scientific); solverInfo.precision(9);
	solverInfo << "Iterations = " << iterations << endl;
	solverInfo << "Residual   = " << residual << endl << endl;

	for (auto const &value : weights)
	{
		solverInfo << value << endl;
	}
	solverInfo.close();
}

//...........................................................................

void MCG::WriteSolverInfoInConsole(const int iterations, const double residual) const
{
	cout.setf(ios::scientific); cout.precision(9);
	cout << "\t Iterations = " << iterations << endl;
	cout << "\t Residual   = " << residual << endl << endl;
	cout.unsetf(ios::scientific);
}

//...........................................................................
