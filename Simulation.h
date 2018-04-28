#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <fstream>

#include "DynamicsStructs.h"

class Simulation
{
public:
	Simulation(DymOptions &dymOptions_in);
	~Simulation();

	int startSimulation();

	/*
	void startSimulation(
		int seed,
		double tempKelvin,
		double impactFactorAu,
		int simulationType
	);
	*/

	void additionalOptions(std::string flag, bool option);

private:
	DymOptions dymOptions_;

	std::ofstream posVel_;

	bool stopSimulation;

	/*
	double timeStep;
	int iterationLoop;
	int printLoop;
	double initialDistance;
	double impactParameter;
	double initialSpeed;
	bool checkStopSimulationConditions;
	bool printEnergy;
	bool simmetrize;
	bool printMovie;
	bool printPosVel;
	double maxStopSimulationDistance;
	double temperatureUnit;
	double mProton;
	std::string outputName;
	*/


	void checkStopSimulation(std::vector<double> &x);

	void printCoulombAtoms(
		std::vector<double> & atoms,
		std::string testName,
		std::vector<double> &atomsCharge);

	void printPositionsAndVelocities(
		std::vector<double> &x,
		std::vector<double> &v,
		std::ofstream &posVelFile_);

//	void initialVelocityKinecticTheory(double TempKelvin);

};

#endif

