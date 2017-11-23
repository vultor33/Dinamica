#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <string>
#include <fstream>

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

	void setOptions(bool printEnergy_in, bool simmetrize_in);

private:
	int integratorType;
	bool printEnergy;
	bool simmetrize;
	std::vector<double> rksParams;
	std::vector<double> atomsCharge;
	std::vector<double> atomsMass;
	
	std::ofstream printEnergyFile_;

};

#endif

