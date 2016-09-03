#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
#include <vector>
#include <cmath>
#include <set>

using namespace std;

class GMRES
{
public:
	vector < set <int> > ig;

	int n, m;
	int m_maxiter, nIter, p;
	vector <int> ia, ja;
	vector <double> di, ggu, ggl;
	vector <double> Mdi, Mggu, Mggl;

	vector <double> weights, f, R0, W, G, C, S;
	vector <vector <double>> H, V;
	double betta, m_eps, oldbetta, cureps;

public:

	GMRES()
		: n(0)
		, m_maxiter(100000)
		, m_eps(1e-09)
	{ };
	~GMRES() { };

	void LinComb(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int n); //�������� �������������� ��������� � Y
	double ScMult(const vector <double> &x, const vector <double> &y, const int n); //��������� ������������
	double NormVect(const vector <double> &x, const int n); //����� �������
	void AVec(const vector <double> &x, const double al, vector <double> &y, const int n); //��������� ������� �� ������
	void dPrintVec(const char *f, const vector <double> &x, const int n); //����� ������� �������������� ����� ������� ��������
	void Ax(const vector <double> &x, vector <double> &b, const int n); //b=A*x            ====��������� ������� �� ������ 
	int Calcx(); //������� ����������� ���� Hy=g
	void Givens(vector <double> &Hi, const int i); //��������� H[*][i] �� ������� ������� � g=G[i]g
	int Solve(); // GMRES
	void AssemblRo(const vector <double> &X, vector <double> &Y); // r0=L(-1)(f-Ax0) === ������ ���
	void ExtractX0(vector <double> &X, vector <double> &Y); // x0=U(-1)r0      === �������� ���
	void MultiplyMatrixOnVector(const vector <double> &x, vector <double> &y) const; // y=U*x           === ������������ ������� �� ������
	void CreateLU();

	// ����� ��. �������
	void PrintMatr(const vector <vector<double>> &A, const int n);
};
