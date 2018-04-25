#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <string>
#include <fstream>

#include "TrajectorySymmetrizer.h"

class Integrator
{
public:
	Integrator();
	~Integrator();

	// nao funciona em campo magnetico
	void rungeKuttaSimetrico(
		std::vector<double> & xInitial,
		std::vector<double> & vInitial,
		double timeStep
		);

	void setAdditionalParams(
		std::vector<double> &atomsMass_in,
		std::vector<double> &atomsCharge_in);

	void setOptions(bool printEnergy_in, int symmetrize_in);

private:
	int integratorType;
	bool printEnergy;
	int symmetrize;
	std::vector<double> rksParams;
	std::vector<double> atomsCharge;
	std::vector<double> atomsMass;
	
	std::ofstream printEnergyFile_;
	TrajectorySymmetrizer symm_;

};

#endif

