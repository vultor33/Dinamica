#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string>
#include <sstream>

#include "Simulation.h"

using namespace std;

void runSimulations(int seedInitial, int seedFinal, double impacFactorAu, double tempKelvin);

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Simlation type not found" << endl;
		cout << "format:  [simultation type (run)]  seedI seedF" << endl;
		exit(1);
	}

	string simlationType = argv[1];
	int seedI, seedF;
	if (simlationType == "run")
	{
		stringstream convert;
		convert << argv[2] << " " << argv[3];
		convert >> seedI >> seedF;
	}

	runSimulations(seedI, seedF, 0.0e0, 100.0e0);
	runSimulations(seedI, seedF, 0.0e0, 300.0e0);
	runSimulations(seedI, seedF, 0.0e0, 1000.0e0);
	runSimulations(seedI, seedF, 0.5e0, 100.0e0);
	runSimulations(seedI, seedF, 0.5e0, 300.0e0);
	runSimulations(seedI, seedF, 0.5e0, 1000.0e0);
	runSimulations(seedI, seedF, 1.0e0, 100.0e0);
	runSimulations(seedI, seedF, 1.0e0, 300.0e0);
	runSimulations(seedI, seedF, 1.0e0, 1000.0e0);
	runSimulations(seedI, seedF, 1.5e0, 100.0e0);
	runSimulations(seedI, seedF, 1.5e0, 300.0e0);
	runSimulations(seedI, seedF, 1.5e0, 1000.0e0);

	/*
	int seed, simulationType;
	double tempKelvin, impactFactorAu;
	seed = time(NULL);
	tempKelvin = 300e0;
	impactFactorAu = 1.0e0;
	simulationType = 0;
	Simulation tacaMagia;
	tacaMagia.startSimulation(seed, tempKelvin, impactFactorAu, simulationType);
	*/
	return 0;
}

void runSimulations(int seedInitial, int seedFinal, double impacFactorAu, double tempKelvin)
{
	for (int i = seedInitial; i <= seedFinal; i++)
	{
		Simulation dim1, dim2, dim3, dim4, dim5;
		dim1.startSimulation(i, tempKelvin, impacFactorAu, 0);
		dim2.startSimulation(i, tempKelvin, impacFactorAu, 1);
		dim3.startSimulation(i, tempKelvin, impacFactorAu, 2);
		dim4.startSimulation(i, tempKelvin, impacFactorAu, 3);
		dim5.startSimulation(i, tempKelvin, impacFactorAu, 4);
	}
}


/*
DIDNT WORKED
MANUALLY REDUCE SPEED
double normV1 = sqrt(v[0] * v[0] + v[4] * v[4] + v[8] * v[8]);
if (normV1 > 10.0e0)
{
v[0] *= 0.99;
v[4] *= 0.99;
v[8] *= 0.99;
}
double normV2 = sqrt(v[2] * v[2] + v[6] * v[6] + v[10] * v[10]);
if (normV2 > 10.0e0)
{
v[2] *= 0.99;
v[6] *= 0.99;
v[10] *= 0.99;
}

BOX
for (int i = 0; i < 12; i++)
{
if ((x[i] > 7.0e0) || (x[i] < -7.0e0))
v[i] *= -1.0e0;
}

*/
