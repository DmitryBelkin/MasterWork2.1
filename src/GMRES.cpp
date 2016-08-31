#include "GMRES.h"

GMRES::GMRES()
{
	int i;

	file1 = fopen("parameters.txt", "rt");
	//считываем размерность пространства
	fscanf(file1, "%d", &n);
	//считываем размерность пространства Крылова
	fscanf(file1, "%d", &m);
	//считываем требуемую точность
	fscanf(file1, "%lg", &eps);
	//считываем максимально допустимое число итераций
	fscanf(file1, "%d", &maxIter);
	fclose(file1);

	//считываем вектор правой части 
	F = new double[n];
	file1 = fopen("vect.txt", "rt");
	for (i = 0; i<n; i++)
		fscanf(file1, "%lg", &F[i]);
	fclose(file1);

	//считываем точное решение
	X3 = new double[n];
	file1 = fopen("linear.txt", "rt");
	for (i = 0; i<n; i++)
		fscanf(file1, "%lg", &X3[i]);
	fclose(file1);

	// ввод матрицы
	file1 = fopen("matr.txt", "rt");

	fscanf(file1, "%d", &n); // размерность матрицы

	//считываем индексный массив ig
	ig = new int[n + 1];
	for (i = 0; i<n + 1; i++)
		fscanf(file1, "%d", &ig[i]);

	int kol = ig[n] - 1;

	//считываем индексный массив jg
	jg = new int[kol];
	for (i = 0; i<kol; i++)
		fscanf(file1, "%d", &jg[i]);

	//считываем нижний треугольник
	ggl = new double[kol];
	for (i = 0; i<kol; i++)
		fscanf(file1, "%lg", &ggl[i]);

	//считываем верхний треугольник
	ggu = new double[kol];
	for (i = 0; i<kol; i++)
		fscanf(file1, "%lg", &ggu[i]);
	//считываем главную диагональ
	di = new double[n];
	for (i = 0; i<n; i++)
		fscanf(file1, "%lg", &di[i]);

	fclose(file1);

	//задаём начальное приближение 
	file1 = fopen("xo.txt", "rt");

	X0 = new double[n];
	for (i = 0; i<n; i++)
		fscanf(file1, "%lg", &X0[i]);

	fclose(file1);
	// выделим память на все оставшееся
	R0 = new double[n];
	W = new double[n];
	G = new double[m + 1];
	C = new double[m];
	S = new double[m];
	H = new double *[m];

	for (i = 0; i<m; i++)
		H[i] = new double[m + 1];

	V = new double *[m + 1];

	for (i = 0; i<m + 1; i++)
		V[i] = new double[n];
}
// дестракт
GMRES :: ~GMRES()
{
	delete[] F;
	delete[] X3;
	delete[] ig;
	delete[] jg;
	delete[] di;
	delete[] ggu;
	delete[] ggl;
	delete[] X0;
	delete[] Mdi;
	delete[] Mggu;
	delete[] Mggl;

	delete[] H;
	delete[] V;
	delete[] G;
	delete[] C;
	delete[] S;
	delete[] R0;
	delete[] W;
}

void GMRES::LUFactor()
{
	int i, j, k;
	Mdi.resize(n);
	Mggl.resize(ig[n] - 1);
	Mggu.resize(ig[n] - 1);

	Mdi[0] = di[0];
	for (i = 1; i<n; i++)
	{
		for (j = ig[i] - 1; j < ig[i + 1] - 1; j++)
		{
			k = jg[j];
			Mggl[j] = ggl[j] - j_k(i, k);
			Mggu[j] = (ggu[j] - j_k(k, i)) / Mdi[k];
		}
		Mdi[i] = di[i] - j_k(i, i);
	}
}

double GMRES::j_k(int j, int k)
{
	int p, q, pj, qj;
	double result = 0.0;
	p = ig[j] - 1;
	q = ig[k] - 1;
	while (p < (ig[j + 1] - 1) && q < (ig[k + 1] - 1))
	{
		pj = jg[p];
		qj = jg[q];
		if (pj < qj)
			p++;
		else
		if (pj > qj)
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

void GMRES::AssemblRo(vector <double> &X, vector <double> &Y)
{
	int i, j, k;
	for (k = 0; k < n; k++)
	{
		Y[k] = 0.0;
		j = ig[k] - 1;
		while (j < ig[k + 1] - 1)
		{
			i = jg[j];
			Y[k] += Mggl[j] * Y[i];
			j++;
		}
		Y[k] = (X[k] - Y[k]) / Mdi[k];
	}
}


void GMRES::ExtractX0(vector <double> &X, vector <double> &Y)
{
	int i, j, l;
	for (i = 0; i<n; i++)
		Y[i] = X[i];
	for (j = n - 1; j > 0; j--)
	{
		l = ig[j + 1];
		while (l > ig[j])
		{
			i = jg[l - 2];
			Y[i] -= Mggu[l - 2] * Y[j];
			l--;
		}
	}
}

void GMRES::LinComb(vector <double> &x, double al, vector <double> &y, vector <double> &rez, int n)
{
	for (int i = 0; i<n; i++)
		rez[i] = y[i] + al*x[i];
}

//скалярное произведение
double GMRES::ScMult(vector <double> &x, vector <double> &y, int n)
{
	double s = 0;
	for (int i = 0; i<n; i++)
		s += x[i] * y[i];
	return s;
}

//норма вектора
double GMRES::NormVect(vector <double> &x, int n)
{
	double s = 0;
	for (int i = 0; i<n; i++)
		s += x[i] * x[i];
	return sqrt(s);
}

void GMRES::Ux(vector <double> &x, vector <double> &y)
{
	int i, j, k;
	for (i = 0; i<n; i++)
		y[i] = x[i];
	for (i = 1; i<n; i++)
	{
		for (j = ig[i] - 1; j < ig[i + 1] - 1; j++)
		{
			k = jg[j];
			y[k] += Mggu[k] * x[i];
		}
	}
}

//умножение вектора на скаляр
void GMRES::AVec(vector <double> &x, double al, vector <double> &y, int n)
{
	for (int i = 0; i<n; i++)
		y[i] = al * x[i];
}

//вывод вектора действительных чисел двойной точности
void GMRES::dPrintVec(char *f, vector <double> &x, int n)
{
	FILE *file = fopen(f, "wt");
	for (int k = 0; k<n; k++)
		fprintf(file, "%12lf\n", x[k]);
	fprintf(file, "\n");
	fclose(file);

}

//умножение матрицы на вектор
void GMRES::Ax(vector <double> &x, vector <double> &b, int n)
{
	int i, j;
	for (i = 0; i<n; i++)
		b[i] = di[i] * x[i];
	for (i = 0; i<n; i++)
	for (j = ig[i] - 1; j<ig[i + 1] - 1; j++)
	{
		b[i] += ggl[j] * x[jg[j]];
		b[jg[j]] += ggu[j] * x[i];
	}
}

//Решение треугольной СЛАУ Hy=g
int GMRES::Calcx()
{
	double s = 1;
	int i, k;
	for (i = 0; i<p; i++)
		s *= H[i][i];
	if (!(s<1e-30))
	{
		G[p - 1] /= H[p - 1][p - 1];
		for (i = p - 2; i >= 0; i--)
		{
			s = 0;
			for (k = i + 1; k<p; k++)
				s += H[k][i] * G[k];
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

void GMRES::PrintMatr(vector <vector <double>> &A, int n)
{
	int i, j;
	for (i = 0; i<n + 1; i++)
	{
		for (j = 0; j<n; j++)
			printf("%lf ", A[j][i]);
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

void Rot(double& x, double& y, double& c, double& s)
{
	double x1;
	x1 = c*x - s*y;
	y = s*x + c*y;
	x = x1;
}

void GMRES::Givens(vector <double> &Hi, int i)
{
	for (int k1 = 0; k1<i; k1++)
		Rot(Hi[k1], Hi[k1 + 1], C[k1], S[k1]);
	CalcSC(Hi[i], Hi[i + 1], C[i], S[i]);
	Rot(G[i], G[i + 1], C[i], S[i]);
}

int GMRES::Solve() //  решатель 
{
	vector <double> Z(n, 0);
	int i, j, k;
	time_t t1, t2;
	file1 = fopen("exit.txt", "w");
	printf("n=%d\n", n);
	fprintf(file1, "n=%d\n", n);
	printf("m=%d\n", m);
	fprintf(file1, "m=%d\n", m);
	printf("Required precision=%le\n", eps);
	fprintf(file1, "Required precision=%le\n", eps);
	printf("maxIter=%d\n", maxIter);
	fprintf(file1, "maxIter=%d\n", maxIter);

	time(&t1); // засекли время входа в цикл

	LUFactor(); // нашли матрицу предобусловливания

	p = m; // выбрали размер подпространства Крылова

	Ax(X0, W, n);				// 
	LinComb(W, (-1), F, Z, n); //  z =(f-Ax0)
	AssemblRo(Z, R0);			//	r=Lz	
	oldbetta = betta = NormVect(R0, n); // 
	ExtractX0(X0, X0);			// x=U(-1)x  
	nIter = 0;					// 

	//BEGIN
	do
	{
		nIter++;				//
		// 
		//  ГМРЕС - способ решать СЛАУ не в нашем пространстве 
		//  а в пространстве меньшей размерности - пространстве Крылова
		//  мы создаем матрицу вращения - она будет у нас обнулять m-n элементов вектора
		//  т.е. поворачивать вектор, и получать вектор с нулевыми решениями(решений) и значениями(пр.части)
		//  когда мы повернем наш вектор, ну и соответсвенно и базис, получим систему меньшей размерности
		//	и решим её, затем опять повернем => получим результат, причем
		//  гораздо быстрее, итеративно и временно, чем не поворачивающие методы - МСГ, ЛОС и т.д.

		for (i = 0; i<m; i++)
		for (j = 0; j<m + 1; j++)
			H[i][j] = 0;			//	матрица поворота - Гивенса

		G[0] = betta;
		for (i = 1; i<m + 1; i++) 	G[i] = 0; // вектор

		AVec(R0, (1 / betta), V[0], n);  // касается поворота 

		for (i = 0; i<m; i++)							// используем пространство Крылова			
		{											// 
			ExtractX0(V[i], W);
			Ax(W, Z, n);
			AssemblRo(Z, W);
			//  поиск матрицы поворота (метод Гивенса)
			for (k = 0; k <= i; k++)
			{
				H[i][k] = ScMult(W, V[k], n);
				LinComb(V[k], (-H[i][k]), W, W, n);
			}//for k

			H[i][i + 1] = NormVect(W, n);

			AVec(W, (1 / H[i][i + 1]), V[i + 1], n);

			Givens(H[i], i);
			if (fabs(G[i + 1])<eps)
			{
				p = i;
				break;
			}
			// 
		}

		if (abs(p) == 0) break;
		// минимизация ||beta*e-h*y||=>y
		// Hy=g ==>y=H(-1)g
		if (Calcx())
		{
			printf("\nIter=%d\n", nIter);
		}

		// xk= xk+ sum(i=1,p; yi*vi)
		// vi - базис подпространства
		// yi - к-ты 
		for (i = 0; i<p; i++)
		{
			ExtractX0(V[i], Z);
			LinComb(Z, G[i], X0, X0, n);
		}

		if (p<m) break; // если наше подпространство оказалось меньше

		// иначе 
		// r0 = L(-1)(f-Axk)
		Ax(X0, W, n);
		LinComb(W, (-1), F, Z, n);
		AssemblRo(Z, R0);
		betta = NormVect(R0, n);
		cureps = betta / oldbetta;
	} while ((cureps > eps) && (nIter < maxIter));	// выход если достигнута какая-то невязка
	// или превышено число итераций
	time(&t2); // засекли время последней операции 
	printf("nIter=%d\n", nIter);
	fprintf(file1, "nIter=%d\n", nIter);
	printf("Real precision=%le\n", cureps);
	fprintf(file1, "Real precision=%le\n", cureps);
	fprintf(file1, "TIME: %ld\n", t2 - t1);
	dPrintVec("res.txt", X0, n);
	for (i = 0; i<n; i++)
		X0[i] = X3[i] - X0[i];
	printf("OtnPogr=%le\n", NormVect(X0, n) / NormVect(X3, n));
	fprintf(file1, "OtnPogr=%le\n", NormVect(X0, n) / NormVect(X3, n));

	if (cureps <= eps)    return  0;
	if (nIter >= maxIter) return -1;
	return -2;
}
