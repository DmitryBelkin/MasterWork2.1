#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
#include <vector>
#include <cmath>
#include <set>

using namespace std;

const string outputPrefix       = "../resources/output/";
const string solverInfoFilename = outputPrefix + "solverInfo.txt";

class GMRES
{
public:
	int m_dimensionOfSLAE;
	int m_gmresDepth;
	int m_currentGmresDepth;
	int m_maxiter;

	vector < set <int> > ig;
	vector <int   > m_ia, m_ja;
	vector <double> m_di , m_ggu , m_ggl ;
	vector <double> m_Mdi, m_Mggu, m_Mggl;

	// @todo ��������� ��� ����������
	vector <double> weights, f, R0, W, G, C, S;
	vector <vector <double>> H, V;
	double betta, m_eps, oldbetta;

public:

	GMRES()
		: m_dimensionOfSLAE(0)
		, m_maxiter(100000)
		, m_eps(1e-09)
	{
		ofstream solverInfo(solverInfoFilename, ios_base::trunc);
	};
	~GMRES() { };

	void   LinearCombination(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int m_dimensionOfSLAE); //�������� �������������� ��������� � Y
	double ScalarProduct    (const vector <double> &x, const vector <double> &y, const int m_dimensionOfSLAE); //��������� ������������
	double NormVector       (const vector <double> &x, const int m_dimensionOfSLAE); //����� �������
	void   MultiplyVectorOnScalar(const vector <double> &x, const double al, vector <double> &y, const int m_dimensionOfSLAE); //��������� ������� �� ������
	void   MultiplyMatrixOnVector(const vector <double> &x, vector <double> &b, const int m_dimensionOfSLAE); //b=A*x            ====��������� ������� �� ������ 
	int    Calcx(); //������� ����������� ���� Hy=g
	void   Givens(vector <double> &Hi, const int i); //��������� H[*][i] �� ������� ������� � g=G[i]g
	void   Solve();
	void   ForwardElimination(const vector <double> &X, vector <double> &Y) const; // r0=L(-1)(f-Ax0) === ������ ���
	void   BackSubstitution  (const vector <double> &X, vector <double> &Y) const; // x0=U(-1)r0      === �������� ���
	void   DecompositionLU();
	void   WriteSolverInfoInFile   (const int iterations, const double residual) const;
	void   WriteSolverInfoInConsole(const int iterations, const double residual) const;
};
