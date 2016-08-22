#include "MCG.h"
#include <assert.h>

//..............................................................................................

void MCG::SetMaxiter(const int maxiter) { m_maxiter = maxiter; }

//..............................................................................................

void MCG::SetEps(const double eps) { m_eps = eps; }

//..............................................................................................

void MCG::Lx(const vector <double> ELf, const vector <double> Df, vector <double> &x, const vector <double> f) const
{
	for (int i = 0; i < n; ++i) x[i] = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = f[i] / Df[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[i] -= ELf[j] * x[ ja[j] ] / Df[i];
	}
}

//..............................................................................................

void MCG::Ux(const vector <double> EUf, vector <double> &x, const vector <double> f) const
{
	for (int i = 0; i < n; ++i) x[i] = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] += f[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[ ja[j] ] -= EUf[j] * x[i];
	}
}

//..............................................................................................

void MCG::LTx(const vector <double> ELf, const vector <double> Df, vector <double> &x, const vector <double> f) const
{
	for (int i = 0; i < n; ++i) x[i] = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = (x[i] + f[i]) / Df[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[ ja[j] ] -= ELf[j] * x[i];
	}
}

//..............................................................................................

void MCG::UTx(const vector <double> EUf, vector <double> &x, const vector <double> f) const
{
	for (int i = 0; i < n; ++i) x[i] = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = f[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
			x[i] -= EUf[j] * x[ ja[j] ];
	}
}

//..............................................................................................

double MCG::ScalarProduct(const vector <double> x, const vector <double> y) const
{
	double temp = 0;
	for (int i = 0; i < n; ++i)
		temp += x[i] * y[i];
	return temp;
}

//..............................................................................................

double MCG::NormVector(const vector <double> x) const {	return sqrt(ScalarProduct(x, x)); }

//..............................................................................................

void MCG::MultMatrixOnVector(const vector <double> EU, const vector <double> EL, const vector <double> D, const vector <double> vect, vector <double> &res) const
{
	for (int i = 0; i < n; ++i) res[i] = 0;
	for (int i = 0; i < n; ++i)
	{
		res[i] = D[i] * vect[i];
		for (int j = ia[i]; j < ia[i + 1]; ++j)
		{
			res[ ja[j] ] += EU[j] * vect[i];
			res[i]       += EL[j] * vect[ja[j]];
		}
	}
}

//..............................................................................................

void MCG::CreateLU()
{
	int size = ia[n], i, j, kj, ki, j1, i0, i1, k;
	double sd, su, sl;
	L.clear();
	U.clear();
	LUdi.clear();

	L.resize(size);
	U.resize(size);
	LUdi.resize(n);
	for (i = 0; i<n; ++i)
	{
		sd = 0;
		i0 = ia[i];
		i1 = ia[i + 1];
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

	cout << "========================> LU was created sucsessfully! <====================== " << endl;
}

//..............................................................................................

void MCG::MCG_LU()
{
	int i, k;
	vector <double> r, z, r_, temp1, temp2;
	double alpha, beta, t, norm;
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
	for (i = 0; i < n; ++i) z[i] = r_[i] = r[i];
	t = NormVector(f);
	norm = NormVector(r) / t;
	cout << "NormVector(r) = " << NormVector(r) << endl;
	cout << "t             = " << t             << endl;
	k = 1;
	while (norm > m_eps && k < m_maxiter)
	{
		Ux(U, temp1, z);
		MultMatrixOnVector(ggu, ggl, di, temp1, temp2);
		Lx(L, LUdi, temp1, temp2);
		LTx(L, LUdi, temp2, temp1);
		MultMatrixOnVector(ggl, ggu, di, temp2, temp1);
		UTx(U, temp2, temp1);
		alpha = ScalarProduct(r, r) / ScalarProduct(temp2, z);
		for (i = 0; i < n; ++i)
		{
			xtch[i] = xtch[i] + alpha*z[i];
			r[i]   -= alpha*temp2[i];
		}
		beta = ScalarProduct(r, r) / ScalarProduct(r_, r_);
		for (i = 0; i < n; ++i)
		{
			z[i]  = r[i] + beta*z[i];
			r_[i] = r[i];
		}
		k++;
		norm = NormVector(r) / t;

		cout << "\r" << norm;
	}
	cout << endl;
	Ux(U, temp1, xtch);
	for (i = 0; i < n; ++i) xtch[i] = temp1[i];

	// @todo ������� ���������� �����
	cout << "k = " << k << endl;
	cout << "norm = " << norm << endl << endl;

	ofstream out("output.txt");
	out << "k = " << k << endl;
	out << "norm = " << norm << endl << endl;
	for(int i = 0; i < xtch.size(); ++i)
	{
		out << xtch[i] << endl;
	}
	out.close();
}
