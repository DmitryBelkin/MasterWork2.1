#include "HelpFunctions.h"
#include <sstream>

//...........................................................................

void GetNumbers(std::vector<int> & result, const std::string & s)
{
	result.clear();
	std::stringstream stream(s);
	int iValue;
	while (stream >> iValue)
		result.push_back(iValue);
}

//...........................................................................
