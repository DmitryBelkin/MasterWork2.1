#include <iostream>
#include <vector>

using namespace std;

class GMRES
{
private:

	int n, m;
	int m_maxiter, nIter, p;
	vector <int> ia, ja;
	// �������� ������� 
	vector <double> di, ggu, ggl;
	// ������� ������������������
	vector <double> Mdi, Mggu, Mggl;

	vector <double> X3, X0, F, R0, W, G, C, S;
	vector <vector <double>> H, V;
	double betta, m_eps, oldbetta, cureps;

public:

	GMRES()
		: n(0)
		, m_maxiter(100000)
		, m_eps(1e-09)
	{ };
	~GMRES() { };

	//�������� �������������� ��������� � Y
	void LinComb(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int n);
	//��������� ������������
	double ScMult(const vector <double> &x, const vector <double> &y, const int n);
	//����� �������
	double NormVect(const vector <double> &x, const int n);

	//��������� ������� �� ������
	void AVec(const vector <double> &x, const double al, vector <double> &y, const int n);
	//����� ������� �������������� ����� ������� ��������
	void dPrintVec(const char *f, const vector <double> &x, const int n);
	//b=A*x            ====��������� ������� �� ������ 
	void Ax(vector <double> &x, vector <double> &b, int n);
	// ����� ��. �������
	void PrintMatr(const vector <vector<double>> &A, const int n);
	//������� ����������� ���� Hy=g
	int Calcx();
	//��������� H[*][i] �� ������� ������� � g=G[i]g
	void Givens(vector <double> &Hi, const int i);
	// GMRES
	int Solve();
	// �������� LU - ����������
	void LUFactor();
	// j-�� ������� �� �-�� ������ -����� ������������ ���������
	double j_k(const int j, const int k);
	// r0=L(-1)(f-Ax0) === ������ ���
	void AssemblRo(const vector <double> &X, vector <double> &Y);
	// x0=U(-1)r0      === �������� ���
	void ExtractX0(vector <double> &X, vector <double> &Y);
	// y=U*x           === ������������ ������� �� ������
	void Ux(vector <double> &x, vector <double> &y);
};
