#include "CalculateScreening.h"

#include <fstream>
#include <iostream>
#include <string>

#include "DynamicsStructs.h"
#include "ReadDymInput.h"
#include "Simulation.h"
#include "Analyze.h"

using namespace std;

CalculateScreening::CalculateScreening(){}

CalculateScreening::~CalculateScreening(){}


void CalculateScreening::screenDynamicCenter(int symmetricType, int initialPositions, int kOnly)
{
	
	int k = 1;
	for (int iEle = -6; iEle <= 8; iEle++)
	{
		for (int iPro = -3; iPro <= 15; iPro++)
		{
			for (int iAngle = 1; iAngle < 90; iAngle++)
			{
				if ((kOnly != -1) && (kOnly != k))
				{
					k++;
					continue;
				}

				double rProton = (double)iPro * 0.1e0 + 0.5156992;
				double rEle = (double)iEle * 0.1e0 + 0.893217217;

				ReadDymInput readDym_;
				readDym_.defineMethodSymmetries(symmetricType, initialPositions, rEle, rProton, (double)iAngle);
				readDym_.addIToName(k);
				k++;
				DymOptions dymOptions_ = readDym_.getDymOptions();
				srand(dymOptions_.seed);

				Simulation sim(dymOptions_);
				int stopStatus = sim.startSimulation();
				if (stopStatus == 1)
				{
					k--;
				}
				else
				{
					Analyze an_;
					an_.chargeDistribution(dymOptions_);
				}

			}
		}
	}

}


void CalculateScreening::screenDynamicRear(int symmetricType, int initialPositions, int kOnly)
{

	int k = 1;
	for (int iEle = -6; iEle <= 9; iEle++)
	{
		for (int iPro = -3; iPro <= 15; iPro++)
		{
			for (int iAngle = 1; iAngle < 90; iAngle++)
			{
				if ((kOnly != -1) && (kOnly != k))
				{
					k++;
					continue;
				}

				double rProton = (double)iPro * 0.1e0 + 0.5156992e0;
				double rEle = (double)iEle * 0.1e0 + rProton + 0.7e0;

				ReadDymInput readDym_;
				readDym_.defineMethodSymmetries(symmetricType, initialPositions, rEle, rProton, (double)iAngle);
				readDym_.addIToName(k);
				k++;
				DymOptions dymOptions_ = readDym_.getDymOptions();
				srand(dymOptions_.seed);

				Simulation sim(dymOptions_);
				int stopStatus = sim.startSimulation();
				if (stopStatus == 1)
				{
					k--;
				}
				else
				{
					Analyze an_;
					an_.chargeDistribution(dymOptions_);
				}

			}
		}
	}

}


void CalculateScreening::screenDynamicCenterPure(int symmetricType, int initialPositions, double iAngle, int kOnly)
{
	int k = 1;
	for (int iEle = -6; iEle <= 8; iEle++)
	{
		for (int iPro = -3; iPro <= 15; iPro++)
		{
			if ((kOnly != -1) && (kOnly != k))
			{
				k++;
				continue;
			}

			double rProton = (double)iPro * 0.1e0 + 0.5156992;
			double rEle = (double)iEle * 0.1e0 + 0.893217217;

			ReadDymInput readDym_;
			readDym_.defineMethodSymmetries(symmetricType, initialPositions, rEle, rProton, iAngle);
			readDym_.addIToName(k);
			k++;
			DymOptions dymOptions_ = readDym_.getDymOptions();
			srand(dymOptions_.seed);

			Simulation sim(dymOptions_);
			int stopStatus = sim.startSimulation();
			if (stopStatus == 1)
			{
				k--;
			}
			else
			{
				Analyze an_;
				an_.chargeDistribution(dymOptions_);
			}
		}
	}

}





void CalculateScreening::screenDynamicRearPure(int symmetricType, int initialPositions, double iAngle, int kOnly)
{
	int k = 1;
	for (int iEle = -6; iEle <= 9; iEle++)
	{
		for (int iPro = -3; iPro <= 15; iPro++)
		{
			if ((kOnly != -1) && (kOnly != k))
			{
				k++;
				continue;
			}

			double rProton = (double)iPro * 0.1e0 + 0.5156992e0;
			double rEle = (double)iEle * 0.1e0 + rProton + 0.7e0;

			ReadDymInput readDym_;
			readDym_.defineMethodSymmetries(symmetricType, initialPositions, rEle, rProton, iAngle);
			readDym_.addIToName(k);
			k++;
			DymOptions dymOptions_ = readDym_.getDymOptions();
			srand(dymOptions_.seed);

			Simulation sim(dymOptions_);
			int stopStatus = sim.startSimulation();
			if (stopStatus == 1)
			{
				k--;
			}
			else
			{
				Analyze an_;
				an_.chargeDistribution(dymOptions_);
			}
		}
	}

}




void  CalculateScreening::calcOne()
{
	int symmetricType = 2;
	int initialPositions = 1;
	int iPro = 0;
	int iEle = 8;
	int iAngle = 89.0e0;
	int k = 0;

	double rProton = (double)iPro * 0.1e0 + 0.5156992;
	double rEle = (double)iEle * 0.1e0 + 0.893217217;

	ReadDymInput readDym_;
	readDym_.defineMethodSymmetries(symmetricType, initialPositions, rEle, rProton, iAngle);
	//readDym_.analyzePlot();
	readDym_.addIToName(k);
	readDym_.electronPlot();
	k++;
	DymOptions dymOptions_ = readDym_.getDymOptions();
	srand(dymOptions_.seed);

	Simulation sim(dymOptions_);
	int stopStatus = sim.startSimulation();
	if (stopStatus == 1)
	{
		k--;
	}
	else
	{
		Analyze an_;
		an_.chargeDistribution(dymOptions_);
	}
}

void CalculateScreening::analyzeAll()
{
	ifstream files_("allFiles.txt");

	string line;
	while (getline(files_, line))
	{
		if (line == "")
			break;

		ReadDymInput readDym_;
		readDym_.defineMethodSymmetries(1, 0, 0.0e0, 0.0e0, 0.0e0);
		DymOptions dymOptions_ = readDym_.getDymOptions();
		dymOptions_.outName = line;
		Analyze an_;
		an_.chargeDistribution(dymOptions_);
	}
}



/* symmetrize
Initial pos = 0
symm - 1 - BH2a (0) -> Ci  -> BH2b  (90)
symm - 2 - LH2a (0) -> C2x -> BH2b  (90)
symm - 3 - BH2a (0) -> C2z -> LH2b' (90)
symm - 4 - LH2a (0) -> Cs  -> LH2b' (90)

Initial pos = 1
symm - 2 - BH2b (0) -> C2x -> LH2b  (90)
*/