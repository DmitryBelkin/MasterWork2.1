#include "FEM.h"

int main()
{
	FEM* problem = new FEM();
	problem->SolveProblem ();
	problem->TransformMeshAfterDisplacement();
	problem->SetBordersOfCheckMesh(-2, 10, -2, 5, -2, 20); // @todo גגמהטעü טח פאיכא
	problem->GenerateMeshForCheck();
	problem->PrintFigure();
	return 0;
}
