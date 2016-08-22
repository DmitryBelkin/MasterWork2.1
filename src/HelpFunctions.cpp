#include "HelpFunctions.h"

void get_numbers(vector<int> & result, const string & s) {
	bool found = false;
	int number = 0;
	bool isMinus = false;
	
	for (string::size_type i = 0; i < s.length(); i++)
	{
		const char ch = s[i];

		if(ch == '-' && s[i+1] >= '0' && s[i+1] <= '9')
		{
			isMinus = true;
		}
		else
		if (ch >= '0' && ch <= '9') 
		{
			const int digit = ch - '0';
			number = number*10 + digit;
			found = true;
		}
		else 
		{
			if (found) 
			{
				if(isMinus == true) { number = number * (-1); }
				result.push_back(number);
				isMinus = false;
				number = 0;
				found = false;
			}
		}
	}
	
	if (found) { result.push_back(number); }
}