#include "CalculateScreening.h"

#include <fstream>
#include <iostream>

#include "DynamicsStructs.h"
#include "ReadDymInput.h"
#include "Simulation.h"
#include "Analyze.h"

using namespace std;

CalculateScreening::CalculateScreening(){}

CalculateScreening::~CalculateScreening(){}


void CalculateScreening::screenDynamicCenter(int symmetricType, int initialPositions, int kOnly)
{
	
	int k = 0;
	for (int iEle = -6; iEle <= 15; iEle++)
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
				bool forcedStop = sim.startSimulation();
				if (forcedStop)
				{
					ofstream excelResult_;
					excelResult_.open(dymOptions_.excelResultsName.c_str(), std::ofstream::out | std::ofstream::app);
					excelResult_ << dymOptions_.outName << ";failed" << endl;
					excelResult_.close();
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

	int k = 0;
	for (int iEle = -6; iEle <= 15; iEle++)
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
				bool forcedStop = sim.startSimulation();
				if (forcedStop)
				{
					ofstream excelResult_;
					excelResult_.open(dymOptions_.excelResultsName.c_str(), std::ofstream::out | std::ofstream::app);
					excelResult_ << dymOptions_.outName << ";failed" << endl;
					excelResult_.close();
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


void  CalculateScreening::calcOne()
{
	int iPro = 0;
	int iEle = 0;
	int iAngle = 90;
	int k = 0;

	double rProton = (double)iPro * 0.1e0 + 0.5156992;
	double rEle = (double)iEle * 0.1e0 + 0.893217217;

	ReadDymInput readDym_;
	readDym_.defineMethodSymmetries(3, 0, rEle, rProton, (double)iAngle);
	readDym_.electronPlot();
	readDym_.addIToName(k);
	k++;
	DymOptions dymOptions_ = readDym_.getDymOptions();
	srand(dymOptions_.seed);

	Simulation sim(dymOptions_);
	bool forcedStop = sim.startSimulation();
	if (forcedStop)
	{
		ofstream excelResult_;
		excelResult_.open(dymOptions_.excelResultsName.c_str(), std::ofstream::out | std::ofstream::app);
		excelResult_ << dymOptions_.outName << ";failed" << endl;
		excelResult_.close();
	}
	else
	{
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