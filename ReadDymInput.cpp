#include "ReadDymInput.h"

#include <string>
#include <sstream>
#include <time.h>

#include "DynamicsStructs.h"

using namespace std;

ReadDymInput::ReadDymInput()
{
	generateDefaultOptions();
}

ReadDymInput::~ReadDymInput(){}


void ReadDymInput::defineMethodSymmetries(
	int symm,
	int initialPositionType,
	double rElec,
	double rProton,
	double angle)
{
	dymOptions_.simulationType = 10;
	dymOptions_.symmetrize = symm;
	dymOptions_.initialPositionType = initialPositionType;
	dymOptions_.initialDistance = 0.0e0;
	dymOptions_.impactParameter = 0.0e0;
	dymOptions_.initialSpeed = 0.0e0;
	dymOptions_.rElec = rElec;
	dymOptions_.rProton = rProton;
	dymOptions_.angle = angle;
}


/*
void ReadDymInput::defineMethodBohr(double angle)
{
	dymOptions_.simulationType = 6;
	dymOptions_.symmetrize = 1;
	dymOptions_.initialDistance = 0.0e0;
	dymOptions_.impactParameter = 0.0e0;
	dymOptions_.initialSpeed = 0.0e0;
	dymOptions_.angleBohrModel = angle;
}

void ReadDymInput::defineMethodEllipseBohr(double rElec, double rProton)
{
	dymOptions_.simulationType = 7;
	dymOptions_.symmetrize = 1;
	dymOptions_.initialDistance = 0.0e0;
	dymOptions_.impactParameter = 0.0e0;
	dymOptions_.initialSpeed = 0.0e0;
	dymOptions_.rElecEllipse = rElec;
	dymOptions_.rProtonEllipse = rProton;
}

void ReadDymInput::defineMethodEllipseBohrAngle(double rElec, double rProton, double angle)
{
	dymOptions_.simulationType = 8;
	dymOptions_.symmetrize = 0;
	dymOptions_.initialDistance = 0.0e0;
	dymOptions_.impactParameter = 0.0e0;
	dymOptions_.initialSpeed = 0.0e0;
	dymOptions_.rElecEllipse = rElec;
	dymOptions_.rProtonEllipse = rProton;
	dymOptions_.angleEllipse = angle;
}

void ReadDymInput::defineMethodLh2b(double rElec, double rProton)
{
	dymOptions_.simulationType = 9;
	dymOptions_.symmetrize = 2;
	dymOptions_.initialDistance = 0.0e0;
	dymOptions_.impactParameter = 0.0e0;
	dymOptions_.initialSpeed = 0.0e0;
	dymOptions_.rElecEllipse = rElec;
	dymOptions_.rProtonEllipse = rProton;
}

*/


void ReadDymInput::electronPlot()
{
	dymOptions_.plotElectronTraject = true;
	dymOptions_.timeStep /= 100.0e0;
}

void ReadDymInput::analyzePlot()
{
	dymOptions_.plotAnalyzeGraphs = true;
	dymOptions_.printPosVel = true;
	dymOptions_.printEnergy = true;
}

void ReadDymInput::setSeed(int seed)
{
	if (seed == -1)
		dymOptions_.seed = (int)time(NULL);
	else
		dymOptions_.seed = seed;
}

void ReadDymInput::addIToName(int i)
{
	stringstream convert;
	convert
		<< i
		<< "-simtype-" << dymOptions_.simulationType
		<< "-initialtype-" << dymOptions_.initialPositionType
		<< "-symmetrizeType-" << dymOptions_.symmetrize
		<< "-rElec-" << dymOptions_.rElec
		<< "-rProton-" << dymOptions_.rProton
		<< "-angle-" << dymOptions_.angle;
	dymOptions_.outName = "simulation-" + convert.str() + ".xyz";
}


DymOptions ReadDymInput::getDymOptions()
{
	return dymOptions_;
}


void ReadDymInput::generateDefaultOptions()
{
	DymOptions dymOptionsDefault;

	dymOptionsDefault.outName = "simulation.xyz";

	dymOptionsDefault.iterationLoop = 10000;
	dymOptionsDefault.printLoop = 300;
	dymOptionsDefault.timeStep = 0.001e0;

	//stop when particcles are too far
	dymOptionsDefault.checkStopSimulationConditions = true;
	dymOptionsDefault.maxStopSimulationDistance = 1.0e6;

	//activate symmetrization
	dymOptionsDefault.symmetrize = 0;

	//print options
	dymOptionsDefault.printEnergy = false;
	dymOptionsDefault.printPosVel = false;
	dymOptionsDefault.printMovie = true;

	//analyze options
	dymOptionsDefault.excelResultsName = "simulationResults.csv";
	dymOptionsDefault.plotAnalyzeGraphs = false;
	dymOptionsDefault.plotElectronTraject = false;

	//initial conditions options
	dymOptionsDefault.seed = 3;
	dymOptionsDefault.simulationType = 0; //simulationType
	dymOptionsDefault.initialDistance = 5.0e0;
	dymOptionsDefault.impactParameter = 0.5e0;
	dymOptionsDefault.initialSpeed = initialVelocityKinecticTheory(300.0e0);
	dymOptionsDefault.energy = -1.17444904341371e0;
	dymOptionsDefault.rElec = 1.0e0;
	dymOptionsDefault.rProton = 0.5e0;
	dymOptionsDefault.angle = 30.0e0;

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



