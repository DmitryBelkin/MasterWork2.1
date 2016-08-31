#include <iostream>
#include <vector>

using namespace std;

class GMRES
{
private:

	int n, m;
	int maxIter, nIter, p;
	vector <int> ig, jg;
	// �������� ������� 
	vector <double> di, ggu, ggl;
	// ������� ������������������
	vector <double> Mdi, Mggu, Mggl;

	vector <double> X3, X0, F, R0, W, G, C, S;
	vector <vector <double>> H, V;
	double betta, eps, oldbetta, cureps;

public:

	GMRES();
	~GMRES();

	//�������� �������������� ��������� � Y
	void LinComb(vector <double> &x, double al, vector <double> &y, vector <double> &rez, int n);
	//��������� ������������
	double ScMult(vector <double> &x, vector <double> &y, int n);
	//����� �������
	double NormVect(vector <double> &x, int n);

	//��������� ������� �� ������
	void AVec(vector <double> &x, double al, vector <double> &y, int n);
	//����� ������� �������������� ����� ������� ��������
	void dPrintVec(vector <char> &f, vector <double> &x, int n);
	//b=A*x            ====��������� ������� �� ������ 
	void Ax(vector <double> &x, vector <double> &b, int n);
	// ����� ��. �������
	void PrintMatr(vector <vector<double>> &A, int n);
	//������� ����������� ���� Hy=g
	int Calcx();
	//��������� H[*][i] �� ������� ������� � g=G[i]g
	void Givens(vector <double> &Hi, int i);
	// GMRES
	int Solve();
	// �������� LU - ����������
	void LUFactor();
	// j-�� ������� �� �-�� ������ -����� ������������ ���������
	double j_k(int j, int k);
	// r0=L(-1)(f-Ax0) === ������ ���
	void AssemblRo(vector <double> &X, vector <double> &Y);
	// x0=U(-1)r0      === �������� ���
	void ExtractX0(vector <double> &X, vector <double> &Y);
	// y=U*x           === ������������ ������� �� ������
	void Ux(vector <double> &x, vector <double> &y);
};
