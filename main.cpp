#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string>
#include <sstream>

#include "DynamicsStructs.h"
#include "Simulation.h"
#include "Analyze.h"

using namespace std;

//void runSimulations(int seedInitial, int seedFinal, double impacFactorAu, double tempKelvin);

//void runSimulationsSymmetric(int seedInitial, int seedFinal, double impacFactorAu, double tempKelvin);

//void generateSeeds(); // GENERATE A LOT OF INITIAL CONDITIONS

DymOptions generateDefaultOptions();

double initialVelocityKinecticTheory(double TempKelvin);

double mProton = 1836.15273443449e0;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		int seed, simulationType;
		double tempKelvin, impactFactorAu;
		/*
		seed = time(NULL);
		tempKelvin = 3000e0;
		impactFactorAu = 0.5e0;
		simulationType = 5;
		*/
		seed = (int)time(NULL);
		tempKelvin = 0.0e0;
		impactFactorAu = 0.0e0;
		simulationType = 6;

		DymOptions dymOptions_ = generateDefaultOptions();
		
		Simulation tacaMagia(dymOptions_);
		tacaMagia.startSimulation(seed, tempKelvin, impactFactorAu, simulationType);

		Analyze an_;
		an_.chargeDistribution("simulacao.xyz");
	}
	else
	{
		/*
		string simlationType = argv[1];
		int seedI, seedF;
		if (simlationType == "run")
		{
			stringstream convert;
			convert << argv[2] << " " << argv[3];
			convert >> seedI >> seedF;
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
		}
		else if (simlationType == "symmetric")
		{
			stringstream convert;
			convert << argv[2] << " " << argv[3];
			convert >> seedI >> seedF;
			runSimulationsSymmetric(seedI, seedF, 0.0e0, 100.0e0);
			runSimulationsSymmetric(seedI, seedF, 0.0e0, 300.0e0);
			runSimulationsSymmetric(seedI, seedF, 0.0e0, 1000.0e0);
			runSimulationsSymmetric(seedI, seedF, 0.5e0, 100.0e0);
			runSimulationsSymmetric(seedI, seedF, 0.5e0, 300.0e0);
			runSimulationsSymmetric(seedI, seedF, 0.5e0, 1000.0e0);
			runSimulationsSymmetric(seedI, seedF, 1.0e0, 100.0e0);
			runSimulationsSymmetric(seedI, seedF, 1.0e0, 300.0e0);
			runSimulationsSymmetric(seedI, seedF, 1.0e0, 1000.0e0);
			runSimulationsSymmetric(seedI, seedF, 1.5e0, 100.0e0);
			runSimulationsSymmetric(seedI, seedF, 1.5e0, 300.0e0);
			runSimulationsSymmetric(seedI, seedF, 1.5e0, 1000.0e0);
		}
		else if (simlationType == "charge")
		{
			Analyze an_;
			an_.takeChargeDistribution();
		}
		*/
	}
	return 0;
}


DymOptions generateDefaultOptions()
{
	DymOptions dymOptionsDefault;

	dymOptionsDefault.outName = "simulacao.xyz";

	//simulation type and steps
	dymOptionsDefault.iterationLoop = 1000;
	dymOptionsDefault.printLoop = 300;
	dymOptionsDefault.timeStep = 0.01e0;

	//stop when particcles is too far
	dymOptionsDefault.checkStopSimulationConditions = true;
	dymOptionsDefault.maxStopSimulationDistance = 100.0e0;

	//activate symmetrization
	dymOptionsDefault.symmetrize = false;

	//print options
	dymOptionsDefault.printEnergy = false;
	dymOptionsDefault.printPosVel = false;
	dymOptionsDefault.printMovie = true;

	//initial conditions options
	dymOptionsDefault.seed = 3;
	dymOptionsDefault.integratorOption = 0; //simulationType
	dymOptionsDefault.initialDistance = 5.0e0;
	dymOptionsDefault.impactParameter = 0.5e0;
	dymOptionsDefault.initialSpeed = initialVelocityKinecticTheory(300.0e0);

	return dymOptionsDefault;

}



double initialVelocityKinecticTheory(double TempKelvin)
{
	// v = sqrt( 3 k t / m) - gas kinetic theory
	double temperatureUnit = 315774.64e0;
	double TempAUnits = TempKelvin / temperatureUnit;
	return sqrt(3.0e0 * TempAUnits / (1.0e0 + mProton));
}


/*

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

void runSimulationsSymmetric(int seedInitial, int seedFinal, double impacFactorAu, double tempKelvin)
{
	for (int i = seedInitial; i <= seedFinal; i++)
	{
		Simulation dim1;
		dim1.startSimulation(i, tempKelvin, impacFactorAu, 5);
	}
}

void generateSeeds()
{
	int nSeeds = 60;
	int nProc = 8;
	ofstream roda_("roda.x");
	roda_ << "#!/bin/bash" << endl << endl;
	for (int i = 480; i < 480 + nSeeds * nProc; i += nSeeds)
	{
		roda_ << "./drkai.x symmetric " << i + 1 << "  " << i + nSeeds << " & " << endl;
	}
	roda_.close();
	system("chmod u+x roda.x");
}

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
