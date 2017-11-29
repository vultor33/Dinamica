#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;

int main()
{
	int nSeeds = 60;
	int nProc = 8;
	ofstream roda_("roda.x");
	roda_ << "#!/bin/bash" << endl << endl;
	for(int i = 480; i < 480 + nSeeds * nProc; i += nSeeds)
	{
		roda_ << "./drkai.x symmetric " << i + 1 << "  " << i + nSeeds << " & " << endl;			
	}
	roda_.close();
	system("chmod u+x roda.x");

	return 0;
}


