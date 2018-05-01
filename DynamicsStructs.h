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
	double maxStopSimulationDistance;
	int symmetrize;
	double toleranceForFinalEnergy;

	//print options
	bool checkStopSimulationConditions;
	bool printEnergy;
	bool printMovie;
	bool printPosVel;

	//analyze options
	std::string excelResultsName;
	bool plotAnalyzeGraphs;
	bool plotElectronTraject;

	//initial conditions options
	int simulationType;
	double initialDistance;
	double impactParameter;
	double initialSpeed;
	double tempKelvin;
	double impactFactorAu;
	double energy;
	double rElec;
	double rProton;
	double angle;
	int initialPositionType; // 0 center ; 1 rear


	//paramters optoins
	double temperatureUnit;
	double mProton;


};


#endif
