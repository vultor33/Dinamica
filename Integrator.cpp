#include "Integrator.h"

#include <vector>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

#include "Fitness.h"

using namespace std;

Integrator::Integrator()
{
	printEnergy = false;
	simmetrize = false;
	integratorType = 0;
	rksParams.resize(19);
	rksParams[0] = 0.095176255;
	rksParams[1] = 0.666296894;
	rksParams[2] = -0.127950286;
	rksParams[3] = 0.024618901;
	rksParams[4] = 0.105972953;
	rksParams[5] = -0.410725534;
	rksParams[6] = 0.448222277;
	rksParams[7] = 0.657729262;
	rksParams[8] = -0.021421199;
	rksParams[9] = -0.875839047;
	rksParams[10] = -0.021421199;
	rksParams[11] = 0.657729262;
	rksParams[12] = 0.448222277;
	rksParams[13] = -0.410725534;
	rksParams[14] = 0.105972953;
	rksParams[15] = 0.024618901;
	rksParams[16] = -0.127950286;
	rksParams[17] = 0.666296894;
	rksParams[18] = 0.095176255;

}

Integrator::~Integrator(){}

void Integrator::setAdditionalParams(
	vector<double> &atomsMass_in,
	vector<double> &atomsCharge_in)
{
	atomsMass = atomsMass_in;
	atomsCharge = atomsCharge_in;
}

void Integrator::setOptions(bool printEnergy_in, bool simmetrize_in)
{
	printEnergy = printEnergy_in;
	if (printEnergy)
	{
		printEnergyFile_.open("printEnergyFile.csv");
	}
	simmetrize = simmetrize_in;
}

void Integrator::rungeKuttaSimetrico(
	std::vector<double> & xInitial,
	std::vector<double> & vInitial,
	double timeStep
)
{
	Fitness fit_;

	if (printEnergy)
	{
		fit_.calculateTotalCoulombSystemEnergy(
			xInitial,
			vInitial,
			atomsCharge,
			atomsMass,
			printEnergyFile_);
		fit_.printCenterOfMass(xInitial, atomsMass, printEnergyFile_);
		printEnergyFile_ << endl;
	}

	size_t size = xInitial.size();

	vector<double> force(size);

	for (size_t i = 0; i < size; i++)
		xInitial[i] += timeStep * rksParams[0] * vInitial[i];


	for (size_t k = 1; k < rksParams.size(); k += 2)
	{
		//fit_.lennardJonesGradient(xInitial, force);
		fit_.CoulombGradient(xInitial, force, atomsCharge);
		for (size_t i = 0; i < size; i++)
		{
			vInitial[i] += timeStep * rksParams[k] * force[i] / atomsMass[i];
			xInitial[i] += timeStep * rksParams[k+1] * vInitial[i];
		}
	}

	if (simmetrize)
	{
		xInitial[2] = -xInitial[0];
		xInitial[3] = -xInitial[1];
		xInitial[6] = -xInitial[4];
		xInitial[7] = -xInitial[5];
		xInitial[10] = -xInitial[8];
		xInitial[11] = -xInitial[9];
		vInitial[2] = -vInitial[0];
		vInitial[3] = -vInitial[1];
		vInitial[6] = -vInitial[4];
		vInitial[7] = -vInitial[5];
		vInitial[10] = -vInitial[8];
		vInitial[11] = -vInitial[9];
	}

}





