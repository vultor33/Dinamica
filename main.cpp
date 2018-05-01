#include "DynamicsStructs.h"
#include "ReadDymInput.h"
#include "Simulation.h"
#include "Analyze.h"
#include "CalculateScreening.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		CalculateScreening calScreen_;

		/*
		Trajetorias com angulo zero removidas, fazer elas depois.
		*/


		//calScreen_.analyzeAll();



		//RUNNING
		//calScreen_.screenDynamicCenter(1, 0); // a-BH2a-BH2b
		//calScreen_.screenDynamicCenter(2, 0); // b-LH2a-BH2b
		//calScreen_.screenDynamicCenter(3, 0); // c-BH2a-LH2bl
		//calScreen_.screenDynamicCenter(4, 0); // d-LH2a-LH2bl
		calScreen_.screenDynamicRear(2, 1); //   e-BH2a-LH2b
		//calScreen_.screenDynamicCenterPure(1, 0, 0.0e0);  // f-BH2a
		//calScreen_.screenDynamicCenterPure(1, 0, 90.0e0);  // g-BH2b - center
		//calScreen_.screenDynamicRearPure(1, 1, 90.0e0); // h-BH2b - rear
		//calScreen_.screenDynamicCenterPure(2, 0, 0.0e0);  // i-LH2a
		//calScreen_.screenDynamicCenterPure(3, 0, 90.0e0);  // j-LH2bl
		//calScreen_.screenDynamicRearPure(2, 1, 0.0e0); // k-LH2b

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
