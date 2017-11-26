#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>
#include <fstream>

class Simulation
{
public:
	Simulation();
	~Simulation();

	void startSimulation(
		int seed,
		double tempKelvin,
		double impactFactorAu,
		int simulationType
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
	bool printMovie;
	bool printPosVel;
	double maxStopSimulationDistance;
	double temperatureUnit;
	double mProton;
	std::string outputName;
	std::ofstream posVel_;

	void checkStopSimulation(std::vector<double> &x);

	void printCoulombAtoms(
		std::vector<double> & atoms,
		std::string testName,
		std::vector<double> &atomsCharge);

	void printPositionsAndVelocities(
		std::vector<double> &x,
		std::vector<double> &v,
		std::ofstream &posVelFile_);

	void initialVelocityKinecticTheory(double TempKelvin);

};

#endif