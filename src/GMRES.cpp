#include "GMRES.h"
#include <numeric>

//...........................................................................

void GMRES::DecompositionLU()
{
	int size = m_ia[m_dimensionOfSLAE], j, kj, ki, j1, k;
	unsigned int i;
	double su, sl;
	m_Mggl.clear();
	m_Mggu.clear();
	m_Mdi.clear();

	m_Mggl.resize(size);
	m_Mggu.resize(size);
	m_Mdi.resize(m_dimensionOfSLAE);
	for (i = 0; i < m_dimensionOfSLAE; ++i)
	{
		double sd = 0;
		const int i0 = m_ia[i];
		const int i1 = m_ia[i + 1];
		for (k = i0; k < i1; ++k)
		{
			su = 0;
			sl = 0;
			j = m_ja[k];
			kj = m_ia[j];
			ki = m_ia[i];
			j1 = m_ia[j + 1];
			while ((ki<k) && (kj<j1))
			{
				if (m_ja[kj] == m_ja[ki])
				{
					sl += m_Mggl[ki] * m_Mggu[kj];
					su += m_Mggl[kj] * m_Mggu[ki];
					ki++; kj++;
				}
				else
				{
					if (m_ja[ki] < m_ja[kj]) ki++;
					else                 kj++;
				}
			}
			assert(m_Mdi[j] != 0);
			m_Mggl[k] =  m_ggl[k] - sl;
			m_Mggu[k] = (m_ggu[k] - su) / m_Mdi[j];
			sd += m_Mggl[k] * m_Mggu[k];
		}
		m_Mdi[i] = m_di[i] - sd;
	}

	cout << "LU - factorization passed" << endl;
}

//...........................................................................

void GMRES::ForwardElimination(const vector <double> &X, vector <double> &Y) const
{
	for (int k = 0; k < m_dimensionOfSLAE; ++k)
	{
		Y[k] = 0.0;
		int j = m_ia[k];
		while (j < m_ia[k + 1])
		{
			const int i = m_ja[j];
			Y[k] += m_Mggl[j] * Y[i];
			j++;
		}
		Y[k] = (X[k] - Y[k]) / m_Mdi[k];
	}
}

//...........................................................................

void GMRES::BackSubstitution(const vector <double> &X, vector <double> &Y) const
{
	for (int i = 0; i < m_dimensionOfSLAE; ++i)
	{
		Y[i] = X[i];
	}
	for (int j = m_dimensionOfSLAE - 1; j > 0; j--)
	{
		int l = m_ia[j + 1];
		while (l > m_ia[j])
		{
			const int i = m_ja[l - 1];
			Y[i] -= m_Mggu[l - 1] * Y[j];
			l--;
		}
	}
}

//...........................................................................

void GMRES::LinearCombination(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int m_dimensionOfSLAE)
{
	for (int i = 0; i < m_dimensionOfSLAE; ++i)
	{
		rez[i] = y[i] + al*x[i];
	}
}

//...........................................................................

//скалярное произведение
double GMRES::ScalarProduct(const vector <double> &x, const vector <double> &y, const int m_dimensionOfSLAE)
{
	return inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

//...........................................................................

//норма вектора
double GMRES::NormVector(const vector <double> &x, const int m_dimensionOfSLAE)
{
	double s = 0;
	for (int i = 0; i < m_dimensionOfSLAE; ++i)
		s += x[i] * x[i];
	return sqrt(s);
}

//...........................................................................

void GMRES::MultiplyVectorOnScalar(const vector <double> &x, const double al, vector <double> &y, const int m_dimensionOfSLAE)
{
	for (int i = 0; i < m_dimensionOfSLAE; ++i)
		y[i] = al * x[i];
}

//...........................................................................

void GMRES::MultiplyMatrixOnVector(const vector <double> &x, vector <double> &b, const int m_dimensionOfSLAE)
{
	for (int i = 0; i < m_dimensionOfSLAE; ++i)
	{
		b[i] = m_di[i] * x[i];
	}
	for (int i = 0; i < m_dimensionOfSLAE; ++i)
	{
		for (int j = m_ia[i]; j < m_ia[i + 1]; ++j) //
		{
			b[ i     ] += m_ggl[j] * x[ m_ja[j] ];
			b[ m_ja[j] ] += m_ggu[j] * x[ i     ];
		}
	}
}

//...........................................................................

//Решение треугольной СЛАУ Hy=g
int GMRES::Calcx()
{
	double s = 1;
	for (int i = 0; i < m_currentGmresDepth; ++i)
	{
		s *= H[i][i];
	}
	if (!(s < 1e-30))
	{
		G[m_currentGmresDepth - 1] /= H[m_currentGmresDepth - 1][m_currentGmresDepth - 1];
		for (int i = m_currentGmresDepth - 2; i >= 0; i--)
		{
			s = 0;
			for (int k = i + 1; k < m_currentGmresDepth; ++k)
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

//...........................................................................

void CalcSC(double& x, double& y, double& c, double& s)
{
	c = x / sqrt(x*x + y*y);
	s = (-y)*c / x;
	x = sqrt(x*x + y*y);
	y = 0;
}

//...........................................................................

void Rot(double& x, double& y, const double& c, const double& s)
{
	const double x1 = c*x - s*y;
	y = s*x + c*y;
	x = x1;
}

//...........................................................................

void GMRES::Givens(vector <double> &Hi, const int i)
{
	for (int k1 = 0; k1 < i; ++k1)
	{
		Rot(Hi[k1], Hi[k1 + 1], C[k1], S[k1]);
	}
	CalcSC(Hi[i], Hi[i + 1], C[i], S[i]);
	Rot(G[i], G[i + 1], C[i], S[i]);
}

//...........................................................................

void GMRES::Solve()
{
	vector <double> Z(m_dimensionOfSLAE, 0);
	DecompositionLU();
	m_currentGmresDepth = m_gmresDepth; // выбрали размер подпространства Крылова
	MultiplyMatrixOnVector(weights, W, m_dimensionOfSLAE);					//
	LinearCombination(W, (-1), f, Z, m_dimensionOfSLAE);			// z =(f-Ax0)
	ForwardElimination(Z, R0);					// r=Lz
	oldbetta = betta = NormVector(R0, m_dimensionOfSLAE); //
	BackSubstitution(weights, weights);		// x=U(-1)x
	int iteration = 0;							// 

	double residual = 0;
	do
	{
		iteration++;
		for (int i = 0; i < m_gmresDepth; ++i)
		{
			for (int j = 0; j < m_gmresDepth + 1; ++j)
			{
				H[i][j] = 0; // матрица поворота - Гивенса
			}
		}

		G[0] = betta;
		for (int i = 1; i < m_gmresDepth + 1; ++i) G[i] = 0;

		MultiplyVectorOnScalar(R0, (1 / betta), V[0], m_dimensionOfSLAE);  // касается поворота 
		for (int i = 0; i < m_gmresDepth; ++i) // используем пространство Крылова
		{
			BackSubstitution(V[i], W);
			MultiplyMatrixOnVector(W, Z, m_dimensionOfSLAE);
			ForwardElimination(Z, W);
			// поиск матрицы поворота (метод Гивенса)
			for (int k = 0; k <= i; ++k)
			{
				H[i][k] = ScalarProduct(W, V[k], m_dimensionOfSLAE);
				LinearCombination(V[k], (-H[i][k]), W, W, m_dimensionOfSLAE);
			}

			H[i][i + 1] = NormVector(W, m_dimensionOfSLAE);
			MultiplyVectorOnScalar(W, (1 / H[i][i + 1]), V[i + 1], m_dimensionOfSLAE);
			Givens(H[i], i);
			if (fabs(G[i + 1]) < m_eps)
			{
				m_currentGmresDepth = i;
				break;
			}
		}

		if (abs(m_currentGmresDepth) == 0) break;
		// минимизация ||beta*e-h*y||=>y
		// Hy=g ==>y=H(-1)g
		if (Calcx())
		{
			printf("\iteration=%d\m_dimensionOfSLAE", iteration);
		}

		// xk= xk+ sum(i=1,m_currentGmresDepth; yi*vi)
		// vi - базис подпространства
		// yi - к-ты 
		for (int i = 0; i < m_currentGmresDepth; ++i)
		{
			BackSubstitution(V[i], Z);
			LinearCombination(Z, G[i], weights, weights, m_dimensionOfSLAE);
		}

		if (m_currentGmresDepth < m_gmresDepth) break; // если наше подпространство оказалось меньше

		// иначе 
		// r0 = L(-1)(f-Axk)
		MultiplyMatrixOnVector(weights, W, m_dimensionOfSLAE);
		LinearCombination(W, (-1), f, Z, m_dimensionOfSLAE);
		ForwardElimination(Z, R0);
		betta = NormVector(R0, m_dimensionOfSLAE);
		residual = betta / oldbetta;
		cout << "\r" << "\t Current residual = " << residual;
	} while ((residual > m_eps) && (iteration < m_maxiter));

	cout << endl << "Solution obtained" << endl;
	if (residual > m_eps)
		cout << "Warning! residual > eps! : residual = " << residual << " eps = " << m_eps << endl;
	if (iteration >= m_maxiter)
		cout << "Warning! Solution obtained by maximum iterations: " << iteration << endl;

	WriteSolverInfoInFile   (iteration, residual);
	WriteSolverInfoInConsole(iteration, residual);
}

//...........................................................................

void GMRES::WriteSolverInfoInFile(const int iterations, const double residual) const
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

void GMRES::WriteSolverInfoInConsole(const int iterations, const double residual) const
{
	cout.setf(ios::scientific); cout.precision(9);
	cout << "\t Iterations = " << iterations << endl;
	cout << "\t Residual   = " << residual << endl << endl;
	cout.unsetf(ios::scientific);
}

//...........................................................................
