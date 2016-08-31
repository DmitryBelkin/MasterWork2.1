#include <iostream>
#include <vector>

using namespace std;

class GMRES
{
private:

	int n, m;
	int m_maxiter, nIter, p;
	vector <int> ia, ja;
	// исходна€ матрица 
	vector <double> di, ggu, ggl;
	// матрица предобусловливани€
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

	//линейное комбинирование результат в Y
	void LinComb(const vector <double> &x, const double al, const vector <double> &y, vector <double> &rez, const int n);
	//скал€рное произведение
	double ScMult(const vector <double> &x, const vector <double> &y, const int n);
	//норма вектора
	double NormVect(const vector <double> &x, const int n);

	//умножение вектора на скал€р
	void AVec(const vector <double> &x, const double al, vector <double> &y, const int n);
	//вывод вектора действительных чисел двойной точности
	void dPrintVec(const char *f, const vector <double> &x, const int n);
	//b=A*x            ====умножение матрицы на вектор 
	void Ax(vector <double> &x, vector <double> &b, int n);
	// вывод кв. матрицы
	void PrintMatr(const vector <vector<double>> &A, const int n);
	//–ешение треугольной —Ћј” Hy=g
	int Calcx();
	//”множение H[*][i] на матрицы √ивенса и g=G[i]g
	void Givens(vector <double> &Hi, const int i);
	// GMRES
	int Solve();
	// неполное LU - разложение
	void LUFactor();
	// j-ый столбец на к-ую строку -сумма произведений элементов
	double j_k(const int j, const int k);
	// r0=L(-1)(f-Ax0) === пр€мой ход
	void AssemblRo(const vector <double> &X, vector <double> &Y);
	// x0=U(-1)r0      === обратный ход
	void ExtractX0(vector <double> &X, vector <double> &Y);
	// y=U*x           === произведение матрицы на вектор
	void Ux(vector <double> &x, vector <double> &y);
};
