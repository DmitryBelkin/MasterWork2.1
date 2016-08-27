#include "FEM.h"
#include "HelpFunctions.h"
#include <fstream>
#include <cassert>
#include <cmath>

using namespace std;

#define DOF             3
#define DOF_ELEM        8
#define LOCAL_DIMENSION ((DOF)*(DOF_ELEM))
const string inputPrefix              = "../resources/input/";
const string inputMeshName            = inputPrefix + "mesh_name.txt";
const string inputSlaeParameters      = inputPrefix + "slae_parameters.txt";
const string inputElastityParameters  = inputPrefix + "elastity_parameters.txt";
const string inputCheckMeshParameters = inputPrefix + "check_mesh_parameters.txt";

static void ReadStuff(const int amountOfLines, ifstream &file)
{
	string str;
	for (int i = 0; i < amountOfLines; ++i) { getline(file, str); }
}

//...........................................................................

void FEM::Input()
{
	InputMesh               ();
	InputSlaeParameters     ();
	InputElastityParameters ();
	InputCheckMeshParameters();
}

//...........................................................................

void FEM::InputMesh(){
	char meshWay[100], meshName[100];

	ifstream mesh_name(inputMeshName);
	mesh_name >> meshName;
	mesh_name.close();

	strcpy_s(meshWay , inputPrefix.c_str());
	strcat_s(meshName, ".unv"             );
	strcat_s(meshWay , meshName           );

	int buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8;
	int num;

	ifstream WorkingArea(meshWay);
	ReadStuff(20, WorkingArea);

	string str;
	vector <int> numbers;
	vector <double> xyzBuf(DOF);
	// считывание координат вершин
	do
	{
		WorkingArea >> xyzBuf[0] >> xyzBuf[1] >> xyzBuf[2];
		m_xyz.push_back(xyzBuf);
		ReadStuff(1, WorkingArea);
		getline(WorkingArea, str); // считал строку с информацией о следующей вершине
		GetNumbers(numbers, str);
		buf1 = numbers[0];
	}
	while(buf1 != -1);
	
	ReadStuff(2, WorkingArea);
	
	// считывание первой строки с информацией о ребре №1
	getline(WorkingArea, str);
	GetNumbers(numbers, str);
	num = numbers[numbers.size() - 1]; // количество считываемых элементов (прямая, треугольник, тетраэдр)
	while(num == 2) // прохожу ребра
	{
		ReadStuff(2, WorkingArea);
		getline(WorkingArea, str);
		GetNumbers(numbers, str);
		num = numbers[numbers.size() - 1]; // количество считываемых элементов (прямая, треугольник, тетраэдр)
	}
	
	while(num == 4) // прохожу прямоугольники
	{
		ReadStuff(1, WorkingArea);
		getline(WorkingArea, str);
		GetNumbers(numbers, str);
		num = numbers[numbers.size()-1]; // количество считываемых элементов (прямая, треугольник, тетраэдр)
	}
	
	vector <int> nvtrBuf(DOF_ELEM);
	while(num == 8) // считываю конечные элементы (параллелипипеды)
	{
		WorkingArea >> buf1 >> buf2 >> buf3 >> buf4 >> buf5 >> buf6 >> buf7 >> buf8;
		nvtrBuf[0] = buf2 - 1;
		nvtrBuf[1] = buf6 - 1;
		nvtrBuf[2] = buf3 - 1;
		nvtrBuf[3] = buf7 - 1;
		nvtrBuf[4] = buf1 - 1;
		nvtrBuf[5] = buf5 - 1;
		nvtrBuf[6] = buf4 - 1;
		nvtrBuf[7] = buf8 - 1;

		m_nvtr.push_back(nvtrBuf);
		ReadStuff(1, WorkingArea);
		getline(WorkingArea, str);
		GetNumbers(numbers, str);
		num = numbers[numbers.size() - 1]; // количество считываемых элементов (прямая, треугольник, тетраэдр)
	}

	ReadStuff(2, WorkingArea);

	while ( ! WorkingArea.eof() )
	{
		size_t founded;
		int amOfEl;
		WorkingArea >> num;
		if (num == -1) // если считал все группы
			break;     // то выйти
		WorkingArea >> buf1 >> buf1 >> buf1 >> buf1 >> buf1 >> buf1 >> buf2;
		amOfEl = buf2;
		ReadStuff(1, WorkingArea);
		getline(WorkingArea, str);

		founded = str.find("nvk1");
		if(founded!=std::string::npos)
		{
			for(int i = 0; i < amOfEl; ++i)
			{
				WorkingArea >> buf1 >> buf2 >> buf1 >> buf1;
				m_nvk1.push_back(buf2 - 1);
			}
		}
		else
		{
			founded=str.find("nvk2_1");
			if(founded!=std::string::npos)
			{
				for(int i = 0; i < amOfEl; ++i)
				{
					WorkingArea >> buf1 >> buf2 >> buf1 >> buf1;
					m_nvk2_1.push_back(buf2 - 1);
				}
			}
			else
			{
				founded=str.find("nvk2_2");
				if(founded!=std::string::npos)
				{
					for(int i = 0; i < amOfEl; ++i)
					{
						WorkingArea >> buf1 >> buf2 >> buf1 >> buf1;
						m_nvk2_2.push_back(buf2 - 1);
					}
				}
			}
		}
	}
	WorkingArea.close();

	cout << "1. Mesh was readed" << endl;
	cout << "\t Mesh info: " << endl;
	cout << "\t Nodes: " << m_xyz.size() << endl;
	cout << "\t Nodes in 1st boundary conditions: " << m_nvk1.size() << endl;
	cout << "\t Nodes in 2d boundary conditions: " << m_nvk2.size() << endl;
	cout << "\t Finite elements: " << m_nvtr.size() << endl;
}

//...........................................................................

void FEM::InputSlaeParameters()
{
	ifstream in(inputSlaeParameters);
	in >> m_eps >> m_maxiter;
	in.close();
}

//...........................................................................

void FEM::InputElastityParameters()
{
	ifstream in(inputElastityParameters);
	in >> m_nu >> m_E;
	in.close();
}

//...........................................................................

void FEM::InputCheckMeshParameters()
{
	ifstream in(inputCheckMeshParameters);
	in >> m_xLeftCheckMesh >> m_xRightCheckMesh >> m_amountOfStepsX;
	in >> m_yLeftCheckMesh >> m_yRightCheckMesh >> m_amountOfStepsY;
	in >> m_zLeftCheckMesh >> m_zRightCheckMesh >> m_amountOfStepsZ;
	in.close();
}

//...........................................................................

void FEM::SolveProblem()
{
	Input();
	SetDefault();
	GenerateMatrixProfle();
	CreateGlobalMatrixAndRightPart();
	MCG_LU();
}

//...........................................................................

void FEM::GenerateMatrixProfle()
{
	// для связей степеней свободы в каждом блоке K (для каждой базисной функции в узле)
	vector <int> degrees(DOF);
	for (unsigned int i = 0; i < m_xyz.size(); ++i)
	{
		degrees[0] = i * DOF;
		degrees[1] = i * DOF + 1;
		degrees[2] = i * DOF + 2;
		AddDegree(degrees);
	}
	for (unsigned int i = 0; i < m_elemAmount; ++i)
	{
		AddNvtr(m_nvtr[i]);
	}
	for (unsigned int i = 1; i < n + 1; ++i)
	{
		ia[i] = ia[i - 1] + ig[i - 1].size();
	}
	ja.resize(ia[n]);
	int i = 0;
	for (auto const &value : ig)
	{
		for (unsigned int k = 0; k < value.size(); ++k, ++i)
		{
			ja[i] = *std::next(value.begin(), k);
		}
	}
	ggl.resize(ia[n]);
	ggu.resize(ia[n]);
	cout << "2. Matrix profle was created" << endl;
}

//...........................................................................

void FEM::AddNvtr(const vector <int> &mtr)
{
	for (int i = 0; i < DOF_ELEM; ++i)
	{
		for (int j = 0; j < DOF_ELEM; ++j)
		{
			for (int shift = 0; shift < DOF; ++shift)
			{
				const int ki = mtr[i] * DOF + shift;
				const int kj = mtr[j] * DOF + shift;
				if (ki > kj)
				{
					ig[ki].insert(kj);
				}
			}
		}
	}
}

//...........................................................................

void FEM::AddDegree(const vector <int> &mtr)
{
	for (int i = 0; i < DOF; ++i)
	{
		for (int j = 0; j < DOF; ++j)
		{
			const int ki = mtr[i];
			const int kj = mtr[j];
			if (ki > kj)
			{
				ig[ki].insert(kj);
			}
		}
	}
}

//...........................................................................

double FEM::X(int k, double x, int numNvtr) const
{
	const double xLeft  = m_xyz[ m_nvtr[numNvtr][0] ][0];
	const double xRight = m_xyz[ m_nvtr[numNvtr][1] ][0];
	const double length = xRight - xLeft;
	double localBasisFunctionValue = 0;
	assert(xRight > xLeft);
	switch(k)
	{
		case 1  : localBasisFunctionValue = (xRight - x) / length; break;
		case 2  : localBasisFunctionValue = (x - xLeft)  / length; break;
		default : assert(false);                                   break;
	}
	return localBasisFunctionValue;
}

//...........................................................................

double FEM::Y(int k, double y, int numNvtr) const
{
	const double yLeft  = m_xyz[ m_nvtr[numNvtr][0] ][1];
	const double yRight = m_xyz[ m_nvtr[numNvtr][2] ][1];
	const double length = yRight - yLeft;
	double localBasisFunctionValue = 0;
	assert(yRight > yLeft);
	switch(k)
	{
		case 1  : localBasisFunctionValue = (yRight - y) / length; break;
		case 2  : localBasisFunctionValue = (y - yLeft)  / length; break;
		default : assert(false);                                   break;
	}
	return localBasisFunctionValue;
}

//...........................................................................

double FEM::Z(int k, double z, int numNvtr) const
{
	const double zLeft  = m_xyz[ m_nvtr[numNvtr][0] ][2];
	const double zRight = m_xyz[ m_nvtr[numNvtr][5] ][2];
	const double length = zRight - zLeft;
	double localBasisFunctionValue = 0;
	assert(zRight > zLeft);
	switch(k)
	{
		case 1  : localBasisFunctionValue = (zRight - z) / length; break;
		case 2  : localBasisFunctionValue = (z - zLeft)  / length; break;
		default : assert(false);                                   break;
	}
	return localBasisFunctionValue;
}

//...........................................................................

double FEM::DerivativeX(int num, int numNvtr) const
{
	const int    k      = ((num - 1) % 2) + 1;
	const double xLeft  = m_xyz[ m_nvtr[numNvtr][0] ][0];
	const double xRight = m_xyz[ m_nvtr[numNvtr][1] ][0];
	const double length = xRight - xLeft;
	double localBasisFunctionDerivative = 0;
	assert(xRight > xLeft);
	switch(k)
	{
		case 1  : localBasisFunctionDerivative = -1 / length; break;
		case 2  : localBasisFunctionDerivative =  1 / length; break;
		default : assert("Error in DerivativeX!");            break;
	}
	return localBasisFunctionDerivative;
}

//...........................................................................

double FEM::DerivativeY(int num, int numNvtr) const
{
	const int    k      = abs((((num - 1) / 2) % 2) + 1);
	const double yLeft  = m_xyz[ m_nvtr[numNvtr][0] ][1];
	const double yRight = m_xyz[ m_nvtr[numNvtr][2] ][1];
	const double length = yRight - yLeft;
	double localBasisFunctionDerivative = 0;
	assert(yRight > yLeft);
	switch(k)
	{
		case 1  : localBasisFunctionDerivative = -1 / length; break;
		case 2  : localBasisFunctionDerivative =  1 / length; break;
		default : assert("Error in DerivativeY!");            break;
	}
	return localBasisFunctionDerivative;
}

//...........................................................................

double FEM::DerivativeZ(int num, int numNvtr) const
{
	const int    k      = abs((num - 1) / 4 + 1);
	const double zLeft  = m_xyz[ m_nvtr[numNvtr][0] ][2];
	const double zRight = m_xyz[ m_nvtr[numNvtr][5] ][2];
	const double length = zRight - zLeft;
	double localBasisFunctionDerivative = 0;
	assert(zRight > zLeft);
	switch(k)
	{
		case 1  : localBasisFunctionDerivative = -1 / length; break;
		case 2  : localBasisFunctionDerivative =  1 / length; break;
		default : assert("Error in DerivativeZ!");            break;
	}
	return localBasisFunctionDerivative;
}

//...........................................................................

double FEM::ThreeLinearFunctionValue(int k, double x, double y, double z, int numNvtr) const
{
	const int a = ((k - 1) % 2) + 1;
	const int b = abs((((k - 1) / 2) % 2) + 1);
	const int c = abs((k - 1) / 4 + 1);
	return X(a, x, numNvtr) * Y(b, y, numNvtr) * Z(c, z, numNvtr);
}

//...........................................................................

double FEM::ThreeLinearFunctionDerivative(int k, int numNvtr) const
{
	const int a = ((k - 1) % 2) + 1;
	const int b = abs((((k - 1) / 2) % 2) + 1);
	const int c = abs((k - 1) / 4 + 1);
	return DerivativeX(a, numNvtr) * DerivativeY(b, numNvtr) * DerivativeZ(c, numNvtr);
}

//...........................................................................

void FEM::CreateGlobalMatrixAndRightPart()
{
	vector < vector <double> > KLocal(LOCAL_DIMENSION);
	vector          <double>   b     (LOCAL_DIMENSION);

	for (int i = 0; i < LOCAL_DIMENSION; ++i)
	{
		KLocal[i].resize(LOCAL_DIMENSION);
	}

	for (unsigned int i = 0; i < m_nvtr.size(); ++i)
	{
		GenerateLocalStiffnessMatrix(i, KLocal);
		AddLocalToGlobal(m_nvtr[i], KLocal, b);
	}

	cout << "3. Matrix and right part without boundary conditions created" << endl;

	BoundaryConditions();

	cout << "4. Boundary conditions considered" << endl;
}

//...........................................................................

void FEM::SetDefault()
{
	m_elemAmount = m_nvtr.size();
	// на каждую вершину приходится по 3 (DOF = 3) смещения (по трём направлениям)
	n = DOF * m_xyz.size();

	ia.clear  ();
	ja.clear  ();
	ig.clear  ();
	f.clear   ();
	di.clear  ();
	weights.clear();
	r.clear   ();
	z.clear   ();
	t.clear   ();
	temp.clear();
	ggl.clear ();
	ggu.clear ();

	ig.resize  (n);
	f.resize   (n);
	di.resize  (n);
	weights.resize(n);
	r.resize   (n);
	z.resize   (n);
	t.resize   (n);
	temp.resize(n);
	ia.resize  (n + 1);

	ia[0] = 0;
}

//...........................................................................

void FEM::GenerateLocalStiffnessMatrix(const int numNvtr, vector < vector <double> > &KLocal) const
{
	vector< vector <double> > LocalBlockK(DOF);
	for (int i = 0; i < DOF; ++i) { LocalBlockK[i].resize(DOF); }

	for(int i = 0; i < DOF_ELEM; ++i)
	{
		for(int j = i; j < DOF_ELEM; ++j)
		{
			GenerateLocalBlockK(i + 1, j + 1, numNvtr, LocalBlockK);
			AddBlockToLocalStiffnessMatrix(i, j, KLocal, LocalBlockK);
		}
	}

	for(int i = 0; i < LOCAL_DIMENSION; ++i)
	{
		for(int j = 0; j < i; ++j)
		{
			KLocal[i][j] = KLocal[j][i];
		}
	}
}

//...........................................................................

void FEM::GenerateLocalBlockK(const int numi, const int numj, const int numNvtr, vector <vector <double> > &LocalBlockK) const
{
	const double dx1 = DerivativeX(numi, numNvtr);
	const double dy1 = DerivativeY(numi, numNvtr);
	const double dz1 = DerivativeZ(numi, numNvtr);

	const double dx2 = DerivativeX(numj, numNvtr);
	const double dy2 = DerivativeY(numj, numNvtr);
	const double dz2 = DerivativeZ(numj, numNvtr);

	//@todo перепроверить
	const double koeff = m_E * (1. - m_nu) / ( (1. + m_nu) * (1. - m_nu) );
	const double a     = m_nu / (1. - m_nu);
	const double b     = (1. - 2. * m_nu) / (2. * (1. - m_nu));

	LocalBlockK[0][0] = dx1 * dx2 + b * (dz1 * dz2 + dy1 * dy2);
	LocalBlockK[1][0] = a * dx1 * dy2 + b * dx2 * dy1          ; LocalBlockK[1][1] = dy1 * dy2 + b * (dz1 * dz2 + dx1 * dx2);
	LocalBlockK[2][0] = a * dx1 * dz2 + b * dx2 * dz1          ; LocalBlockK[2][1] = a * dy1 * dz2 + b * dy2 * dz1          ; LocalBlockK[2][2] = dz1 * dz2 + b * (dx1 * dx2 + dy1 * dy2);

	// @todo проверить интеграл. В данном случае из-за линыйности бф - считается просто объём элемента
	for (int i = 0; i < DOF; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			LocalBlockK[i][j] *= VolumeOfParallelepiped(numNvtr) * koeff;
			LocalBlockK[j][i] = LocalBlockK[i][j];
		}
	}
}

//...........................................................................

double FEM::VolumeOfParallelepiped(const int numNvtr) const
{
	const double xLength = m_xyz[ m_nvtr[numNvtr][1] ][0] - m_xyz[ m_nvtr[numNvtr][0] ][0];
	const double yLength = m_xyz[ m_nvtr[numNvtr][2] ][1] - m_xyz[ m_nvtr[numNvtr][0] ][1];
	const double zLength = m_xyz[ m_nvtr[numNvtr][4] ][2] - m_xyz[ m_nvtr[numNvtr][0] ][2];

	return xLength * yLength * zLength;
}

//...........................................................................

void FEM::AddBlockToLocalStiffnessMatrix(const int numi, const int numj, vector <vector <double> > &KLocal, const vector <vector <double> > &LocalBlockK) const
{
	for (int i = 0; i < DOF; ++i)
	{
		for (int j = i; j < DOF; ++j)
		{
			assert( ((numi * DOF + i) < LOCAL_DIMENSION) && ((numj * DOF + j) < LOCAL_DIMENSION) );
			KLocal[numi * DOF + i][numj * DOF + j] = LocalBlockK[i][j];
		}
	}
}

//...........................................................................

void FEM::AddLocalToGlobal(const vector <int> &mtrx, const vector < vector <double> > &K, const vector <double> &b)
{
	int k;
	for (int i = 0; i < DOF_ELEM; ++i)
	{
		for (int shift = 0; shift < DOF; ++shift)
		{
			const int ki = mtrx[i] * DOF + shift;

			assert(ki < di.size()); di[ki] += K[i * DOF + shift][i * DOF + shift];
			assert(ki < f.size()) ;  f[ki] += b[i * DOF + shift]                 ;
			for (int j = 0; j < i; ++j)
			{
				int kj = mtrx[j] * DOF + shift;
				if (ki > kj) { assert(ki < ia.size()); k = ia[ki]; }
				else
				{
					assert(kj < ia.size()); 
					k  = ia[kj];
					kj = ki    ;
				}

				while (ja[k] != kj)
				{
					++k;
					assert(k < ja.size());
				}
				assert(k < ggl.size());
				ggl[k] += K[i * DOF + shift][j * DOF + shift];
				ggu[k] += K[i * DOF + shift][j * DOF + shift];
			}
		}
	}
}

//...........................................................................

void FEM::BoundaryConditions()
{
	int uzel;
	double Ug;

	// вторые краевые на одной грани (сверху)
	for (unsigned int i = 0; i < m_nvk2_1.size(); ++i)
	{
		uzel = m_nvk2_1[i];
		Ug = GetTraction_1(m_xyz[uzel][0], m_xyz[uzel][1], m_xyz[uzel][2]);
		f[uzel * DOF + 2] += Ug;
	}

	// вторые краевые на одной грани (снизу)
	for (unsigned int i = 0; i < m_nvk2_2.size(); ++i)
	{
		uzel = m_nvk2_2[i];
		Ug = GetTraction_2(m_xyz[uzel][0], m_xyz[uzel][1], m_xyz[uzel][2]);
		f[uzel * DOF + 2] += Ug;
	}

	for (unsigned int i = 0; i < m_nvk1.size(); ++i)
	{
		uzel = m_nvk1[i];
		Ug = GetUg(m_xyz[uzel][0], m_xyz[uzel][1], m_xyz[uzel][2]);
		for (int shift = 0; shift < DOF; ++shift)
		{
			Boundary_1(uzel * DOF + shift, Ug);
		}
	}
}

//...........................................................................

void FEM::Boundary_1(int num, double Ug)
{
	int j = 0;
	di[num] = 1.0;
	f [num] = Ug;
	for (; j < ia[n]; ++j)
	{
		if (ja[j] == num)
			ggu[j] = 0;
	}
	for (j = ia[num]; j < ia[num + 1]; ++j)
		ggl[j] = 0;
}

//...........................................................................

// величина смещения
double FEM::GetUg(double x, double y, double z)
{
	x; y; z;
	return 0.;
}

//...........................................................................

// величина тяги верхней грани
double FEM::GetTraction_1(double x, double y, double z)
{
	x; y; z;
	return 1.;
}

//...........................................................................

// величина тяги нижней грани
double FEM::GetTraction_2(double x, double y, double z)
{
	x; y; z;
	return -1.;
}

//...........................................................................

void FEM::PrintFigure()
{
	ofstream Figure("Figure_tecplot.dat");
	Figure << "TITLE = \"U\"" << endl;
	Figure << "VARIABLES = \"x\", \"y\", \"z\", \"p\"" << endl;
	//Figure << "ZONE i=" << m_amountOfStepsX << 
	//             ", j=" << m_amountOfStepsY << 
	//             ", k=" << m_amountOfStepsZ << 
	//             ", F=POINT" << endl;
	Figure << "ZONE i=" << 2 << 
	             ", j=" << 2 << 
	             ", k=" << 2 << 
	             ", F=POINT" << endl;

	float value = 0;
	vector <double> xPoints, yPoints, zPoints;
	for (unsigned int i = 0; i < m_checkMesh.size(); ++i)
	{
		if(PointBelongsToArea(m_checkMesh[i][0], m_checkMesh[i][1], m_checkMesh[i][2])) 
			value = 1;
		else 
			value = 0;

		if(value == 1)
		{
			xPoints.push_back(m_checkMesh[i][0]);
			yPoints.push_back(m_checkMesh[i][1]);
			zPoints.push_back(m_checkMesh[i][2]);
			/*Figure << m_checkMesh[i][0] << '\t'
			       << m_checkMesh[i][1] << '\t'
			       << m_checkMesh[i][2] << '\t'
			       << value             << endl;*/
		}
	}

	//сортируем
	sort(xPoints.begin(),xPoints.end());
	//и удаляем дубликаты
	xPoints.resize(unique(xPoints.begin(),xPoints.end())-xPoints.begin());

	//сортируем
	sort(yPoints.begin(),yPoints.end());
	//и удаляем дубликаты
	yPoints.resize(unique(yPoints.begin(),yPoints.end())-yPoints.begin());

	//сортируем
	sort(zPoints.begin(),zPoints.end());
	//и удаляем дубликаты
	zPoints.resize(unique(zPoints.begin(),zPoints.end())-zPoints.begin());
	
	Figure << xPoints[0] << '\t'
	       << yPoints[0] << '\t'
	       << zPoints[0] << '\t'
	       << value             << endl;

	Figure << xPoints[xPoints.size() - 1] << '\t'
	       << yPoints[0] << '\t'
	       << zPoints[0] << '\t'
	       << value             << endl;

	Figure << xPoints[0] << '\t'
	       << yPoints[yPoints.size() - 1] << '\t'
	       << zPoints[0] << '\t'
	       << value             << endl;

	Figure << xPoints[xPoints.size() - 1] << '\t'
	       << yPoints[yPoints.size() - 1] << '\t'
	       << zPoints[0] << '\t'
	       << value             << endl;

	// z 

		Figure
			<< xPoints[0]                  << '\t'
			<< yPoints[0]                  << '\t'
			<< zPoints[zPoints.size() - 1] << '\t'
			<< value                       << endl;

	Figure 
			<< xPoints[xPoints.size() - 1] << '\t'
			<< yPoints[0]                  << '\t'
			<< zPoints[zPoints.size() - 1] << '\t'
			<< value                       << endl;

	Figure 
			<< xPoints[0]                  << '\t'
			<< yPoints[yPoints.size() - 1] << '\t'
			<< zPoints[zPoints.size() - 1] << '\t'
			<< value                       << endl;

	Figure
			<< xPoints[xPoints.size() - 1] << '\t'
			<< yPoints[yPoints.size() - 1] << '\t'
			<< zPoints[zPoints.size() - 1] << '\t'
			<< value                       << endl;

	Figure.close();
}

//...........................................................................

void FEM::TransformMeshAfterDisplacement()
{
	for (unsigned int i = 0; i < m_xyz.size(); ++i)
	{
		m_xyz[i][0] = m_xyz[i][0] + weights[i * DOF    ];
		m_xyz[i][1] = m_xyz[i][1] + weights[i * DOF + 1];
		m_xyz[i][2] = m_xyz[i][2] + weights[i * DOF + 2];
	}
}

//...........................................................................

void FEM::GenerateMeshForCheck()
{
	// @todo сделать шаги членами класса?
	const double xStep = (m_xRightCheckMesh - m_xLeftCheckMesh) / static_cast<double>(m_amountOfStepsX);
	const double yStep = (m_yRightCheckMesh - m_yLeftCheckMesh) / static_cast<double>(m_amountOfStepsY);
	const double zStep = (m_zRightCheckMesh - m_zLeftCheckMesh) / static_cast<double>(m_amountOfStepsZ);

	vector <double> bufCheckMesh(DOF);
	for (unsigned int k = 0; k < m_amountOfStepsZ; ++k)
	{
		for (unsigned int j = 0; j < m_amountOfStepsY; ++j)
		{
			for (unsigned int i = 0; i < m_amountOfStepsX; ++i)
			{
				bufCheckMesh[0] = m_xLeftCheckMesh + xStep * i;
				bufCheckMesh[1] = m_yLeftCheckMesh + yStep * j;
				bufCheckMesh[2] = m_zLeftCheckMesh + zStep * k;
				m_checkMesh.push_back(bufCheckMesh);
			}
		}
	}
}

//...........................................................................

bool FEM::PointBelongsToArea(double x, double y, double z) const
{
	bool ifBelongs = false;
	for (unsigned int i = 0; i < m_nvtr.size() && !ifBelongs; ++i)
	{
		if(PointBelongsToParallelepiped(x, y, z, i))
			ifBelongs = true;
	}
	return ifBelongs;
}

//...........................................................................

bool FEM::PointBelongsToParallelepiped(double x, double y, double z, int numNvtr) const
{
	if(    x >= m_xyz[ m_nvtr[numNvtr][0] ][0] && x <= m_xyz[ m_nvtr[numNvtr][1] ][0]
		&& y >= m_xyz[ m_nvtr[numNvtr][0] ][1] && y <= m_xyz[ m_nvtr[numNvtr][2] ][1]
		&& z >= m_xyz[ m_nvtr[numNvtr][0] ][2] && z <= m_xyz[ m_nvtr[numNvtr][4] ][2])
			return true;

	return false;
}

//...........................................................................
