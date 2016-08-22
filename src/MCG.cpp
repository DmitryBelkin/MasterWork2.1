#include "MCG.h"
#include <assert.h>

//..............................................................................................

void MCG::setMaxiter(int maxiter)
{
	m_maxiter = maxiter;
}

//..............................................................................................

void MCG::setEps(double eps)
{
	m_eps = eps;
}

//..............................................................................................

void MCG::Lx(vector <double> &ELf, vector <double> &Df, vector <double> &x, vector <double> &f)
{
	for (int i = 0; i<n; i++) x[i] = 0;
	for (int i = 0; i<n; i++)
	{
		x[i] = f[i] / Df[i];
		for (int j = ia[i] ; j<ia[i + 1] ; j++) x[i] = x[i] - ELf[j] * x[ja[j] ] / Df[i];
	}
}

//..............................................................................................

void MCG::Ux(vector <double> &EUf, vector <double> &x, vector <double> &f)
{
	for (int i = 0; i<n; i++) x[i] = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = x[i] + f[i];
		for (int j = ia[i] ; j<ia[i + 1] ; j++) x[ja[j] ] = x[ja[j] ] - EUf[j] * x[i];
	}
}

//..............................................................................................

void MCG::LTx(vector <double> &ELf, vector <double> &Df, vector <double> &x, vector <double> &f)
{
	for (int i = 0; i<n; i++) x[i] = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = (x[i] + f[i]) / Df[i];
		for (int j = ia[i]; j<ia[i + 1] ; j++) x[ja[j]] = x[ja[j]] - ELf[j] * x[i];
	}
}

//..............................................................................................

void MCG::UTx(vector <double> &EUf, vector <double> &x, vector <double> &f)
{
	for (int i = 0; i<n; i++) x[i] = 0;
	for (int i = 0; i<n; i++)
	{
		x[i] = f[i];
		for (int j = ia[i]; j<ia[i + 1]; j++) x[i] = x[i] - EUf[j] * x[ja[j] ];
	}
}

//..............................................................................................

double MCG::scalarProduct(vector <double> &x, vector <double> &y)
{
	double temp = 0;
	for (int i = 0; i < n; i++) temp += x[i] * y[i];
	return temp;
}

//..............................................................................................

double MCG::normVector(vector <double> &x)
{
	return sqrt(scalarProduct(x, x));
}

//..............................................................................................

void MCG::multMatrixOnVector(vector <double> &EU, vector <double> &EL, vector <double> &D, vector <double> &vect, vector <double> &res)
{
	for (int i = 0; i<n; i++) res[i] = 0;
	for (int i = 0; i<n; i++)
	{
		res[i] = D[i] * vect[i];
		for (int j = ia[i] ; j<ia[i + 1] ; j++)
		{
			res[ja[j] ] = res[ja[j]] + EU[j] * vect[i];
			res[i] = res[i] + EL[j] * vect[ja[j]];
		}
	}
}

//..............................................................................................

void MCG::createLU()
{
	int size = ia[n], i, j, kj, ki, j1, i0, i1, k;
	double sd, su, sl;
	L.clear();
	U.clear();
	LUdi.clear();

	L.resize(size);
	U.resize(size);
	LUdi.resize(n);
	for (i = 0; i<n; i++)
	{
		sd = 0;
		i0 = ia[i];
		i1 = ia[i + 1];
		for (k = i0; k<i1; k++)
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
					if (ja[ki]<ja[kj]) ki++;
					else kj++;
				}
			}
			assert(LUdi[j] != 0);
			L[k] = ggl[k] - sl;
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
	r.resize(n);
	z.resize(n);
	r_.resize(n);
	temp1.resize(n);
	temp2.resize(n);
	createLU();
	Lx(L, LUdi, temp1, f);
	LTx(L, LUdi, temp2, temp1);
	multMatrixOnVector(ggl, ggu, di, temp2, temp1);
	UTx(U, r, temp1);
	for (i = 0; i<n; i++) z[i] = r_[i] = r[i];
	t = normVector(f);
	norm = normVector(r) / t;
	cout << "normVector(r) = " << normVector(r) << endl;
	cout << "t             = " << t             << endl;
	k = 1;
	while (norm>m_eps && k<m_maxiter)
	{
		Ux(U, temp1, z);
		multMatrixOnVector(ggu, ggl, di, temp1, temp2);
		Lx(L, LUdi, temp1, temp2);
		LTx(L, LUdi, temp2, temp1);
		multMatrixOnVector(ggl, ggu, di, temp2, temp1);
		UTx(U, temp2, temp1);
		alpha = scalarProduct(r, r) / scalarProduct(temp2, z);
		for (i = 0; i<n; i++)
		{
			xtch[i] = xtch[i] + alpha*z[i];
			r[i] = r[i] - alpha*temp2[i];
		}
		beta = scalarProduct(r, r) / scalarProduct(r_, r_);
		for (i = 0; i<n; i++)
		{
			z[i] = r[i] + beta*z[i];
			r_[i] = r[i];
		}
		k++;
		norm = normVector(r) / t;

		cout << "\r" << norm;
	}
	cout << endl;
	Ux(U, temp1, xtch);
	for (i = 0; i < n; i++) xtch[i] = temp1[i];

	// @todo сделать нормальный вывод
	cout << "k = " << k << endl;
	cout << "norm = " << norm << endl << endl;

	ofstream out("output.txt");
	out << "k = " << k << endl;
	out << "norm = " << norm << endl << endl;
	for(int i = 0; i < xtch.size(); i++)
	{
		out << xtch[i] << endl;
	}
	out.close();
}
