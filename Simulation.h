#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <fstream>

#include "DynamicsStructs.h"
#include "Fitness.h"

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
	bool stopSimulation;

	void checkStopSimulation(std::vector<double> &x);

	void printCoulombAtoms(
		std::vector<double> & atoms,
		std::string testName,
		std::vector<double> &atomsCharge);

	void printSimulationInfo(
		std::vector<double> &x,
		std::vector<double> &v,
		std::vector<double> &atomsCharge,
		std::vector<double> &atomsMass,
		std::ofstream &posVelFile_);

	DymOptions dymOptions_;
	std::ofstream posVel_;
	Fitness fit_;

};

#endif

