#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>

class Simulation
{
public:
	Simulation();
	~Simulation();

	void startSimulation(
		int seed
	);

private:
	double timeStep;
	int iterationLoop;
	int printLoop;
	double initialDistance;
	double impactParameter;
	double initialSpeed;
	bool checkStopSimulationConditions;
	bool stopSimulation;
	bool printEnergy;
	bool simmetrize;
	std::string outputName;

	void checkStopSimulation(std::vector<double> &x);

	void printCoulombAtoms(
		std::vector<double> & atoms,
		std::string testName,
		std::vector<double> &atomsCharge);


};

#endif