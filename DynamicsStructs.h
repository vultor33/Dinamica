#ifndef DYNAMICSSTRUCTS_H
#define DYNAMICSSTRUCTS_H

#include <string>


// isso tem q ser como se fosse uma leitura de input

struct DymOptions
{
	std::string outName;

	//integrator options
	double timeStep;
	int iterationLoop;
	int printLoop;
	int seed;
	int integratorOption;//simulationType
	double maxStopSimulationDistance;

	//print options
	bool checkStopSimulationConditions;
	bool printEnergy;
	bool symmetrize;
	bool printMovie;
	bool printPosVel;

	//initial conditions options
	double initialDistance;
	double impactParameter;
	double initialSpeed;
	double tempKelvin, impactFactorAu;

	//paramters optoins
	double temperatureUnit;
	double mProton;


};

#endif
