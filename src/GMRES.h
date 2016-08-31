#include <iostream>
#include <vector>

using namespace std;

class GMRES
{
private:

	int n, m;
	int maxIter, nIter, p;
	vector <int> ig, jg;
	// исходна€ матрица 
	vector <double> di, ggu, ggl;
	// матрица предобусловливани€
	vector <double> Mdi, Mggu, Mggl;

	vector <double> X3, X0, F, R0, W, G, C, S;
	vector <vector <double>> H, V;
	double betta, eps, oldbetta, cureps;

public:

	GMRES();
	~GMRES();

	//линейное комбинирование результат в Y
	void LinComb(vector <double> &x, double al, vector <double> &y, vector <double> &rez, int n);
	//скал€рное произведение
	double ScMult(vector <double> &x, vector <double> &y, int n);
	//норма вектора
	double NormVect(vector <double> &x, int n);

	//умножение вектора на скал€р
	void AVec(vector <double> &x, double al, vector <double> &y, int n);
	//вывод вектора действительных чисел двойной точности
	void dPrintVec(vector <char> &f, vector <double> &x, int n);
	//b=A*x            ====умножение матрицы на вектор 
	void Ax(vector <double> &x, vector <double> &b, int n);
	// вывод кв. матрицы
	void PrintMatr(vector <vector<double>> &A, int n);
	//–ешение треугольной —Ћј” Hy=g
	int Calcx();
	//”множение H[*][i] на матрицы √ивенса и g=G[i]g
	void Givens(vector <double> &Hi, int i);
	// GMRES
	int Solve();
	// неполное LU - разложение
	void LUFactor();
	// j-ый столбец на к-ую строку -сумма произведений элементов
	double j_k(int j, int k);
	// r0=L(-1)(f-Ax0) === пр€мой ход
	void AssemblRo(vector <double> &X, vector <double> &Y);
	// x0=U(-1)r0      === обратный ход
	void ExtractX0(vector <double> &X, vector <double> &Y);
	// y=U*x           === произведение матрицы на вектор
	void Ux(vector <double> &x, vector <double> &y);
};
