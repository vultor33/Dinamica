#include "ReadDymInput.h"

#include <string>

#include "DynamicsStructs.h"

using namespace std;

ReadDymInput::ReadDymInput(){}

ReadDymInput::~ReadDymInput(){}

DymOptions ReadDymInput::generateDefaultOptions()
{
	DymOptions dymOptionsDefault;

	dymOptionsDefault.outName = "simulacao.xyz";

	//simulation type and steps
	/*
	dymOptionsDefault.iterationLoop = 1000;
	dymOptionsDefault.printLoop = 300;
	dymOptionsDefault.timeStep = 0.01e0;
	*/
	dymOptionsDefault.iterationLoop = 1000;
	dymOptionsDefault.printLoop = 300;
	dymOptionsDefault.timeStep = 0.01e0;

	//stop when particcles is too far
	dymOptionsDefault.checkStopSimulationConditions = true;
	dymOptionsDefault.maxStopSimulationDistance = 100.0e0;

	//activate symmetrization
	dymOptionsDefault.symmetrize = false;

	//print options
	dymOptionsDefault.printEnergy = true;
	dymOptionsDefault.printPosVel = true;
	dymOptionsDefault.printMovie = true;

	//initial conditions options
	dymOptionsDefault.seed = 3;
	dymOptionsDefault.simulationType = 0; //simulationType
	dymOptionsDefault.initialDistance = 5.0e0;
	dymOptionsDefault.impactParameter = 0.5e0;
	dymOptionsDefault.initialSpeed = initialVelocityKinecticTheory(300.0e0);
	dymOptionsDefault.energy = -1.17444904341371e0;
	dymOptionsDefault.angleBohrModel = 5.0e0;

	return dymOptionsDefault;

}



double ReadDymInput::initialVelocityKinecticTheory(double TempKelvin)
{
	// v = sqrt( 3 k t / m) - gas kinetic theory
	double mProton = 1836.15273443449e0;
	double temperatureUnit = 315774.64e0;
	double TempAUnits = TempKelvin / temperatureUnit;
	return sqrt(3.0e0 * TempAUnits / (1.0e0 + mProton));
}



