#include "GMRES.h"

void GMRES::CreateLU()
{
	int size = ia[n], j, kj, ki, j1, k;
	unsigned int i;
	double su, sl;
	Mggl.clear();
	Mggu.clear();
	Mdi.clear();

	Mggl.resize(size);
	Mggu.resize(size);
	Mdi.resize(n);
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
					sl += Mggl[ki] * Mggu[kj];
					su += Mggl[kj] * Mggu[ki];
					ki++; kj++;
				}
				else
				{
					if (ja[ki] < ja[kj]) ki++;
					else                 kj++;
				}
			}
			assert(Mdi[j] != 0);
			Mggl[k] =  ggl[k] - sl;
			Mggu[k] = (ggu[k] - su) / Mdi[j];
			sd += Mggl[k] * Mggu[k];
		}
		Mdi[i] = di[i] - sd;
	}

	cout << "LU - factorization passed" << endl;
}

void GMRES::LUFactor()
{
	Mdi.resize(n);
	Mggl.resize(ia[n] - 1);
	Mggu.resize(ia[n] - 1);
	Mdi[0] = di[0];
	for (int i = 1; i < n; ++i)
	{
		for (int j = ia[i]; j < ia[i + 1]; ++j)
		{
			const int k = ja[j];
			Mggl[j] =  ggl[j] - j_k(i, k);
			Mggu[j] = (ggu[j] - j_k(k, i)) / Mdi[k];
		}
		Mdi[i] = di[i] - j_k(i, i);
	}
}

double GMRES::j_k(const int j, const int k)
{
	double result = 0.0;
	int p = ia[j];
	int q = ia[k];
	while (p < (ia[j + 1]) && q < (ia[k + 1]))
	{
		const int pj = ja[p];
		const int qj = ja[q];
		if (pj < qj) p++;
		else 
		if (pj > qj) q++;
		else
		{
			result += Mggl[p] * Mggu[q];
			p++;
			q++;
		}
	}
	return result;
}

void GMRES::AssemblRo(const vector <double> &X, vector <double> &Y)
{
	for (int k = 0; k < n; ++k)
	{
		Y[k] = 0.0;
		int j = ia[k];
		while (j < ia[k + 1])
		{
			const int i = ja[j];
			Y[k] += Mggl[j] * Y[i];
			j++;
		}
		Y[k] = (X[k] - Y[k]) / Mdi[k];
	}
}

void GMRES::ExtractX0(vector <double> &X, vector <double> &Y)
{
	for (int i = 0; i < n; ++i)
	{
		Y[i] = X[i];
	}
	for (int j = n - 1; j > 0; j--)
	{
		int l = ia[j + 1];
		while (l > ia[j])
		{
			const int i = ja[l - 1];
			Y[i] -= Mggu[l - 1] * Y[j];
			l--;
		}
	}
}

void GMRES::LinComb(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int n)
{
	for (int i = 0; i < n; ++i)
	{
		rez[i] = y[i] + al*x[i];
	}
}

//скал€рное произведение
double GMRES::ScMult(const vector <double> &x, const vector <double> &y, const int n)
{
	double s = 0;
	for (int i = 0; i < n; ++i)
		s += x[i] * y[i];
	return s;
}

//норма вектора
double GMRES::NormVect(const vector <double> &x, const int n)
{
	double s = 0;
	for (int i = 0; i < n; ++i)
		s += x[i] * x[i];
	return sqrt(s);
}

void GMRES::Ux(vector <double> &x, vector <double> &y)
{
	for (int i = 0; i < n; ++i)
	{
		y[i] = x[i];
	}
	for (int i = 1; i < n; ++i)
	{
		for (int j = ia[i] - 1; j < ia[i + 1] - 1; ++j)
		{
			const int k = ja[j];
			y[k] += Mggu[k] * x[i];
		}
	}
}

//умножение вектора на скал€р
void GMRES::AVec(const vector <double> &x, const double al, vector <double> &y, const int n)
{
	for (int i = 0; i < n; ++i)
		y[i] = al * x[i];
}

//вывод вектора действительных чисел двойной точности
void GMRES::dPrintVec(const  char *f, const  vector <double> &x, const  int n)
{
	FILE *file = fopen(f, "wt");
	for (int k = 0; k < n; ++k)
		fprintf(file, "%12lf\n", x[k]);
	fprintf(file, "\n");
	fclose(file);
}

void GMRES::Ax(const vector <double> &x, vector <double> &b, const int n)
{
	for (int i = 0; i < n; ++i)
	{
		b[i] = di[i] * x[i];
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = ia[i]; j < ia[i + 1]; ++j) //
		{
			b[ i     ] += ggl[j] * x[ ja[j] ];
			b[ ja[j] ] += ggu[j] * x[ i     ];
		}
	}
}

//–ешение треугольной —Ћј” Hy=g
int GMRES::Calcx()
{
	double s = 1;
	for (int i = 0; i < p; ++i)
	{
		s *= H[i][i];
	}
	if (!(s < 1e-30))
	{
		G[p - 1] /= H[p - 1][p - 1];     //
		for (int i = p - 2; i >= 0; i--) //
		{
			s = 0;
			for (int k = i + 1; k < p; ++k)
			{
				s += H[k][i] * G[k];
			}
			G[i] = (G[i] - s) / H[i][i];
		}
		return 0;
	}
	else
	{
		printf("Matrix is incorrect!");
		return 1;
	}
}

void GMRES::PrintMatr(const vector <vector <double>> &A, const int n)
{
	for (int i = 0; i < n + 1; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			printf("%lf ", A[j][i]);
		}
		printf("\n");
	}
}

void CalcSC(double& x, double& y, double& c, double& s)
{
	c = x / sqrt(x*x + y*y);
	s = (-y)*c / x;
	x = sqrt(x*x + y*y);
	y = 0;
}

void Rot(double& x, double& y, const double& c, const double& s)
{
	const double x1 = c*x - s*y;
	y = s*x + c*y;
	x = x1;
}

void GMRES::Givens(vector <double> &Hi, const int i)
{
	for (int k1 = 0; k1 < i; ++k1)
	{
		Rot(Hi[k1], Hi[k1 + 1], C[k1], S[k1]);
	}
	CalcSC(Hi[i], Hi[i + 1], C[i], S[i]);
	Rot(G[i], G[i + 1], C[i], S[i]);
}

int GMRES::Solve() //  решатель 
{
	vector <double> Z(n, 0);

	CreateLU();
	//LUFactor(); // нашли матрицу предобусловливани€
	p = m; // выбрали размер подпространства  рылова

	Ax(weights, W, n);					//
	LinComb(W, (-1), f, Z, n);			// z =(f-Ax0)
	AssemblRo(Z, R0);					// r=Lz
	oldbetta = betta = NormVect(R0, n); //
	ExtractX0(weights, weights);		// x=U(-1)x
	nIter = 0;							// 

	do
	{
		nIter++;

		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < m + 1; ++j)
			{
				H[i][j] = 0; // матрица поворота - √ивенса
			}
		}

		G[0] = betta;
		for (int i = 1; i < m + 1; ++i) G[i] = 0; // вектор

		AVec(R0, (1 / betta), V[0], n);  // касаетс€ поворота 

		for (int i = 0; i < m; ++i) // используем пространство  рылова
		{
			ExtractX0(V[i], W);
			Ax(W, Z, n);
			AssemblRo(Z, W);
			// поиск матрицы поворота (метод √ивенса)
			for (int k = 0; k <= i; ++k)
			{
				H[i][k] = ScMult(W, V[k], n);
				LinComb(V[k], (-H[i][k]), W, W, n);
			}

			H[i][i + 1] = NormVect(W, n);
			AVec(W, (1 / H[i][i + 1]), V[i + 1], n);
			Givens(H[i], i);
			if (fabs(G[i + 1]) < m_eps)
			{
				p = i;
				break;
			}
		}

		if (abs(p) == 0) break;
		// минимизаци€ ||beta*e-h*y||=>y
		// Hy=g ==>y=H(-1)g
		if (Calcx())
		{
			printf("\nIter=%d\n", nIter);
		}

		// xk= xk+ sum(i=1,p; yi*vi)
		// vi - базис подпространства
		// yi - к-ты 
		for (int i = 0; i < p; ++i)
		{
			ExtractX0(V[i], Z);
			LinComb(Z, G[i], weights, weights, n);
		}

		if (p < m) break; // если наше подпространство оказалось меньше

		// иначе 
		// r0 = L(-1)(f-Axk)
		Ax(weights, W, n);
		LinComb(W, (-1), f, Z, n);
		AssemblRo(Z, R0);
		betta = NormVect(R0, n);
		cureps = betta / oldbetta;
	} while ((cureps > m_eps) && (nIter < m_maxiter));	// выход если достигнута кака€-то нев€зка

	// или превышено число итераций
	//printf("nIter=%d\n", nIter);
	//printf("Real precision=%le\n", cureps);
	//dPrintVec("res.txt", weights, n);
	//for (int i = 0; i < n; ++i)
	//{
	//	weights[i] = X3[i] - weights[i];
	//}
	//printf("OtnPogr=%le\n", NormVect(weights, n) / NormVect(X3, n));

	if (cureps <= m_eps    ) return  0;
	if (nIter  >= m_maxiter) return -1;
	return -2;
}
