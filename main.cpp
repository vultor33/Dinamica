#include <time.h>

#include "DynamicsStructs.h"
#include "ReadDymInput.h"
#include "Simulation.h"
#include "Analyze.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		ReadDymInput readDym_;
		DymOptions dymOptions_ = readDym_.generateDefaultOptions();
		dymOptions_.seed = (int)time(NULL);
		srand(dymOptions_.seed);
		dymOptions_.simulationType = 6;
		dymOptions_.timeStep = 0.0001;

		Simulation tacaMagia(dymOptions_);
		tacaMagia.startSimulation();

		Analyze an_;
		an_.chargeDistribution(dymOptions_.outName);
	}
	return 0;
}


/*
gnuplot
cd "C:\\Users\\frederico\\source\\repos\\Dinamica\\Dinamica"
*/






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
