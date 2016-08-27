#include "FEM.h"

int main()
{
	FEM* problem = new FEM();
	problem->SolveProblem();
	problem->TransformMeshAfterDisplacement();
	problem->GenerateMeshForCheck();
	problem->PrintFigure();
	return 0;
}
