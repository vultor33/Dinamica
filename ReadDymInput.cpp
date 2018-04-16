#include "ReadDymInput.h"

#include <string>
#include <sstream>

#include "DynamicsStructs.h"

using namespace std;

ReadDymInput::ReadDymInput()
{
	generateDefaultOptions();
}

ReadDymInput::~ReadDymInput(){}

DymOptions ReadDymInput::getDymOptions()
{
	return dymOptions_;
}

void ReadDymInput::addIToName(int i)
{
	stringstream convert;
	convert << i;
	dymOptions_.outName = "simulation-" + convert.str() + ".xyz";
}


void ReadDymInput::generateDefaultOptions()
{
	DymOptions dymOptionsDefault;

	dymOptionsDefault.outName = "simulation.xyz";

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

	//analyze options
	dymOptionsDefault.excelResultsName = "simulationResults.csv";
	dymOptionsDefault.plotAnalyzeGraphs = false;

	//initial conditions options
	dymOptionsDefault.seed = 3;
	dymOptionsDefault.simulationType = 0; //simulationType
	dymOptionsDefault.initialDistance = 5.0e0;
	dymOptionsDefault.impactParameter = 0.5e0;
	dymOptionsDefault.initialSpeed = initialVelocityKinecticTheory(300.0e0);
	dymOptionsDefault.energy = -1.17444904341371e0;
	dymOptionsDefault.angleBohrModel = 20.0e0;

	dymOptions_ = dymOptionsDefault;

}



double ReadDymInput::initialVelocityKinecticTheory(double TempKelvin)
{
	// v = sqrt( 3 k t / m) - gas kinetic theory
	double mProton = 1836.15273443449e0;
	double temperatureUnit = 315774.64e0;
	double TempAUnits = TempKelvin / temperatureUnit;
	return sqrt(3.0e0 * TempAUnits / (1.0e0 + mProton));
}



