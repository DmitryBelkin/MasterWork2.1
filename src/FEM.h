#include "MCG.h"

class FEM : MCG
{
public:
	FEM()
		: m_elemAmount(0)
		, m_xLeftCheckMesh(0)
		, m_xRightCheckMesh(0)
		, m_yLeftCheckMesh(0)
		, m_yRightCheckMesh(0)
		, m_zLeftCheckMesh(0)
		, m_zRightCheckMesh(0)
		, m_amountOfStepsX(0)
		, m_amountOfStepsY(0)
		, m_amountOfStepsZ(0)
		, m_nu(0)
		, m_E(0)
	{ }
	virtual ~FEM() { };

	// attributes
	unsigned  int                          m_elemAmount;
	std::vector < std::vector < double > > m_xyz, m_xyzDisplacement, m_checkMesh;
	std::vector < std::vector < int > >    m_nvtr;
	std::vector < int >                    m_nvk1, m_nvk2;
	double                                 m_xLeftCheckMesh, m_xRightCheckMesh,
	                                       m_yLeftCheckMesh, m_yRightCheckMesh,
	                                       m_zLeftCheckMesh, m_zRightCheckMesh;
	unsigned int                           m_amountOfStepsX, m_amountOfStepsY, m_amountOfStepsZ;
	std::vector < int >                    m_nvk2_1, m_nvk2_2;
	double                                 m_nu, m_E;
	vector < vector < double > >           m_C;

	// methods
	void   Input();
	void   InputMesh();
	void   InputSlaeParameters();
	void   InputElastityParameters();

	void   SolveProblem();
	void   SetDefault();
	void   GenerateMatrixProfle();
	void   CreateGlobalMatrixAndRightPart();
	void   GenerateLocalStiffnessMatrix(int numNvtr, vector <vector <double> > &KLocal) const;
	void   GenerateLocalBlockK(int i, int j, int numNvtr, vector <vector <double> > &LocalBlockK) const;
	double VolumeOfParallelepiped(int numNvtr) const;
	bool   PointBelongsToArea(double x, double y, double z) const;
	bool   PointBelongsToParallelepiped(double x, double y, double z, int numNvtr) const;
	void   AddLocalToGlobal(vector <int> &mtrx, vector <vector <double> > &K, vector <double> &b);
	void   AddBlockToLocalStiffnessMatrix(int numi, int numj, vector <vector <double> > &KLocal , vector <vector <double> > &LocalBlockK) const;
	void   AddNvtr(vector <int> &mtr);
	void   AddDegree(vector <int> &mtr);
	void   PrintFigure();
	void   BoundaryConditions();
	void   Boundary_1(int num, double Ug);
	double GetUg(double x, double y, double z);
	double GetTraction_1(double x, double y, double z);
	double GetTraction_2(double x, double y, double z);
	double X(int k, double x, int numNvtr) const;
	double Y(int k, double y, int numNvtr) const;
	double Z(int k, double z, int numNvtr) const;
	double DerivativeX(int num, int numNvtr) const;
	double DerivativeY(int num, int numNvtr) const;
	double DerivativeZ(int num, int numNvtr) const;
	double ThreeLinearFunctionValue(int k, double x, double y, double z, int numNvtr) const;
	double ThreeLinearFunctionDerivative(int k, int numNvtr) const;
	void   TransformMeshAfterDisplacement();
	void   GenerateMeshForCheck();
	void   SetBordersOfCheckMesh(double xLeftCheckMesh, double xRightCheckMesh, double yLeftCheckMesh, double yRightCheckMesh, double zLeftCheckMesh, double zRightCheckMesh);
};
