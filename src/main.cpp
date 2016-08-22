#include "FEM.h"
#include <iostream>

int main(int argc, char *argv[])
{
	FEM* problem = new FEM();
	problem->SolveProblem ();
	problem->TransformMeshAfterDisplacement();
	problem->SetBordersOfCheckMesh(-2, 10, -2, 5, -2, 20); // @todo גגמהטעü טח פאיכא
	problem->GenerateMeshForCheck();
	problem->PrintFigure();
}
