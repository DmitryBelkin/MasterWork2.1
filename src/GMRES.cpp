#include "GMRES.h"

////GMRES::GMRES()
////{
////	//��������� ����������� ������������ �������
////	//m = 
////
////	//��������� ������ ������ ����� 
////	F.resize(n, 0);
////
////	//��������� ������ �������
////	X3.resize(n, 0);
////
////	// ���� �������
////
////	//��������� ��������� ������ ig
////	ig.resize(n + 1, 0);
////
////	int kol = ig[n] - 1;
////
////	//��������� ��������� ������ jg
////	jg.resize(kol, 0);
////
////	//��������� ������ �����������
////	ggl.resize(kol, 0);
////
////	//��������� ������� �����������
////	ggu.resize(kol, 0);
////
////	//��������� ������� ���������
////	di.resize(n, 0);
////
////	//����� ��������� ����������� 
////	X0.resize(n, 0);
////
////	// ������� ������ �� ��� ����������
////	R0.resize(n, 0);
////	W.resize(n, 0);
////	G.resize(m + 1, 0);
////	C.resize(m, 0);
////	S.resize(m, 0);
////	H.resize(m);
////	for (int i = 0; i < m; ++i)
////		H[i].resize(m + 1, 0);
////
////	V.resize(m + 1);
////	for (int i = 0; i < m + 1; ++i)
////		V[i].resize(n, 0);
////}

void GMRES::LUFactor()
{
	Mdi.resize(n);
	Mggl.resize(ig[n] - 1);
	Mggu.resize(ig[n] - 1);
	Mdi[0] = di[0];
	for (int i = 1; i < n; ++i)
	{
		for (int j = ig[i] - 1; j < ig[i + 1] - 1; ++j)
		{
			const int k = jg[j];
			Mggl[j] = ggl[j] - j_k(i, k);
			Mggu[j] = (ggu[j] - j_k(k, i)) / Mdi[k];
		}
		Mdi[i] = di[i] - j_k(i, i);
	}
}

double GMRES::j_k(const int j, const int k)
{
	double result = 0.0;
	int p = ig[j] - 1;
	int q = ig[k] - 1;
	while (p < (ig[j + 1] - 1) && q < (ig[k + 1] - 1))
	{
		const int pj = jg[p];
		const int qj = jg[q];
		if (pj < qj)
			p++;
		else if (pj > qj)
				q++;
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
	int j;
	for (int k = 0; k < n; ++k)
	{
		Y[k] = 0.0;
		j = ig[k] - 1;
		while (j < ig[k + 1] - 1)
		{
			const int i = jg[j];
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
		int l = ig[j + 1];
		while (l > ig[j])
		{
			const int i = jg[l - 2];
			Y[i] -= Mggu[l - 2] * Y[j];
			l--;
		}
	}
}

void GMRES::LinComb(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int n)
{
	for (int i = 0; i < n; ++i)
		rez[i] = y[i] + al*x[i];
}

//��������� ������������
double GMRES::ScMult(const vector <double> &x, const vector <double> &y, const int n)
{
	double s = 0;
	for (int i = 0; i < n; ++i)
		s += x[i] * y[i];
	return s;
}

//����� �������
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
		for (int j = ig[i] - 1; j < ig[i + 1] - 1; ++j)
		{
			const int k = jg[j];
			y[k] += Mggu[k] * x[i];
		}
	}
}

//��������� ������� �� ������
void GMRES::AVec(const vector <double> &x, const double al, vector <double> &y, const int n)
{
	for (int i = 0; i<n; ++i)
		y[i] = al * x[i];
}

//����� ������� �������������� ����� ������� ��������
void GMRES::dPrintVec(const  char *f, const  vector <double> &x, const  int n)
{
	FILE *file = fopen(f, "wt");
	for (int k = 0; k < n; ++k)
		fprintf(file, "%12lf\n", x[k]);
	fprintf(file, "\n");
	fclose(file);

}

//��������� ������� �� ������
void GMRES::Ax(vector <double> &x, vector <double> &b, int n)
{
	for (int i = 0; i < n; ++i)
	{
		b[i] = di[i] * x[i];
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = ig[i] - 1; j < ig[i + 1] - 1; ++j)
		{
			b[i] += ggl[j] * x[jg[j]];
			b[jg[j]] += ggu[j] * x[i];
		}
	}
}

//������� ����������� ���� Hy=g
int GMRES::Calcx()
{
	double s = 1;
	for (int i = 0; i < p; ++i)
	{
		s *= H[i][i];
	}
	if (!(s < 1e-30))
	{
		G[p - 1] /= H[p - 1][p - 1];
		for (int i = p - 2; i >= 0; i--)
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

int GMRES::Solve() //  �������� 
{
	vector <double> Z(n, 0);
	time_t t1, t2;

	LUFactor(); // ����� ������� ������������������
	p = m; // ������� ������ ��������������� �������

	Ax(X0, W, n);				// 
	LinComb(W, (-1), F, Z, n); //  z =(f-Ax0)
	AssemblRo(Z, R0);			//	r=Lz	
	oldbetta = betta = NormVect(R0, n); // 
	ExtractX0(X0, X0);			// x=U(-1)x  
	nIter = 0;					// 

	do
	{
		nIter++;
		// 
		//  ����� - ������ ������ ���� �� � ����� ������������ 
		//  � � ������������ ������� ����������� - ������������ �������
		//  �� ������� ������� �������� - ��� ����� � ��� �������� m-n ��������� �������
		//  �.�. ������������ ������, � �������� ������ � �������� ���������(�������) � ����������(��.�����)
		//  ����� �� �������� ��� ������, �� � ������������� � �����, ������� ������� ������� �����������
		//	� ����� �, ����� ����� �������� => ������� ���������, ������
		//  ������� �������, ���������� � ��������, ��� �� �������������� ������ - ���, ��� � �.�.

		for (int i = 0; i < m; ++i)
			for (int j = 0; j < m + 1; ++j)
				H[i][j] = 0; // ������� �������� - �������

		G[0] = betta;
		for (int i = 1; i < m + 1; ++i) G[i] = 0; // ������

		AVec(R0, (1 / betta), V[0], n);  // �������� �������� 

		for (int i = 0; i < m; ++i) // ���������� ������������ �������
		{
			ExtractX0(V[i], W);
			Ax(W, Z, n);
			AssemblRo(Z, W);
			// ����� ������� �������� (����� �������)
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
		// ����������� ||beta*e-h*y||=>y
		// Hy=g ==>y=H(-1)g
		if (Calcx())
		{
			printf("\nIter=%d\n", nIter);
		}

		// xk= xk+ sum(i=1,p; yi*vi)
		// vi - ����� ���������������
		// yi - �-�� 
		for (int i = 0; i < p; ++i)
		{
			ExtractX0(V[i], Z);
			LinComb(Z, G[i], X0, X0, n);
		}

		if (p < m) break; // ���� ���� ��������������� ��������� ������

		// ����� 
		// r0 = L(-1)(f-Axk)
		Ax(X0, W, n);
		LinComb(W, (-1), F, Z, n);
		AssemblRo(Z, R0);
		betta = NormVect(R0, n);
		cureps = betta / oldbetta;
	} while ((cureps > m_eps) && (nIter < m_maxiter));	// ����� ���� ���������� �����-�� �������

	// ��� ��������� ����� ��������
	printf("nIter=%d\n", nIter);
	printf("Real precision=%le\n", cureps);
	dPrintVec("res.txt", X0, n);
	for (int i = 0; i < n; ++i)
	{
		X0[i] = X3[i] - X0[i];
	}
	printf("OtnPogr=%le\n", NormVect(X0, n) / NormVect(X3, n));

	if (cureps <= m_eps    ) return  0;
	if (nIter  >= m_maxiter) return -1;
	return -2;
}
