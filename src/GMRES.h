#include<stdio.h>
#include<math.h>
#include <time.h>
#include <conio.h>


class GMRES
{
private:

	int n;
	int *ig;
	int *jg;
	// �������� ������� 
	double *di;
	double *ggu;
	double *ggl;
	// ������� ������������������
	double *Mdi;
	double *Mggu;
	double *Mggl;

	int m;

	double *X3, *X0, *F, *R0, *W, *G, **H, **V, *C, *S;
	double betta, eps, oldbetta, cureps;
	int maxIter, nIter, p;
	FILE *file1;

public:

	GMRES();
	~GMRES();


	//�������� �������������� ��������� � Y
	void LinComb(double *x, double al, double *y, double *rez, int n);
	//��������� ������������
	double ScMult(double *x, double *y, int n);
	//����� �������
	double NormVect(double *x, int n);
	//��������� ������� �� ������
	void AVec(double *x, double al, double *y, int n);
	//����� ������� �������������� ����� ������� ��������
	void dPrintVec(char *f, double *x, int n);
	//b=A*x            ====��������� ������� �� ������ 
	void Ax(double *x, double *b, int n);
	// ����� ��. �������
	void PrintMatr(double **A, int n);
	//������� ����������� ���� Hy=g
	int Calcx();
	//��������� H[*][i] �� ������� ������� � g=G[i]g
	void Givens(double *Hi, int i);
	// GMRES
	int Solve();
	// �������� LU - ����������
	void LUFactor();
	// j-�� ������� �� �-�� ������ -����� ������������ ���������
	double j_k(int j, int k);
	// r0=L(-1)(f-Ax0) === ������ ���
	void AssemblRo(double *X, double *Y);
	// x0=U(-1)r0      === �������� ���
	void ExtractX0(double *X, double *Y);
	// y=U*x           === ������������ ������� �� ������
	void Ux(double *x, double *y);
};
