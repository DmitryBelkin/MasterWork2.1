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
	// исходна€ матрица 
	double *di;
	double *ggu;
	double *ggl;
	// матрица предобусловливани€
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


	//линейное комбинирование результат в Y
	void LinComb(double *x, double al, double *y, double *rez, int n);
	//скал€рное произведение
	double ScMult(double *x, double *y, int n);
	//норма вектора
	double NormVect(double *x, int n);
	//умножение вектора на скал€р
	void AVec(double *x, double al, double *y, int n);
	//вывод вектора действительных чисел двойной точности
	void dPrintVec(char *f, double *x, int n);
	//b=A*x            ====умножение матрицы на вектор 
	void Ax(double *x, double *b, int n);
	// вывод кв. матрицы
	void PrintMatr(double **A, int n);
	//–ешение треугольной —Ћј” Hy=g
	int Calcx();
	//”множение H[*][i] на матрицы √ивенса и g=G[i]g
	void Givens(double *Hi, int i);
	// GMRES
	int Solve();
	// неполное LU - разложение
	void LUFactor();
	// j-ый столбец на к-ую строку -сумма произведений элементов
	double j_k(int j, int k);
	// r0=L(-1)(f-Ax0) === пр€мой ход
	void AssemblRo(double *X, double *Y);
	// x0=U(-1)r0      === обратный ход
	void ExtractX0(double *X, double *Y);
	// y=U*x           === произведение матрицы на вектор
	void Ux(double *x, double *y);
};
