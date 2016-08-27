#include "HelpFunctions.h"
#include <sstream>

//...........................................................................

void GetNumbers(std::vector<int> & result, const std::string & s)
{
	result.clear();
	
	std::stringstream stream(s);

	int n;
	while (stream >> n)
		result.push_back(n);
}

//...........................................................................
