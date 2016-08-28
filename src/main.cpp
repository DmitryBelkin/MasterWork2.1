#include "FEM.h"
#include <boost\shared_ptr.hpp>

int main()
{
	std::shared_ptr<FEM> pProblem(new FEM());
	pProblem->SolveProblem();
	pProblem->TransformMeshAfterDisplacement();
	//pProblem->GenerateMeshForCheck();
	//pProblem->PrintFigure();
	pProblem->PrintMesh();
	return 0;
}
