#include "FEM.h"
#include <iostream>

int main(int argc, char *argv[])
{

	FEM* problem = new FEM();
	problem->SolveProblem ();

	problem->TransformMeshAfterDisplacement();
	//                              0   7   0  2   0  16
	problem->SetBordersOfCheckMesh(-2, 10, -2, 5, -2, 20); // @todo גגמהטעü טח פאיכא
	problem->GenerateMeshForCheck();
	problem->PrintFigure();
}