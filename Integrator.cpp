#include "Integrator.h"

#include <vector>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <boost/numeric/odeint.hpp>

#include "Fitness.h"
#include "TrajectorySymmetrizer.h"
#include "FitnessCoulombOdeint.h"

typedef std::vector< double > state_type;
typedef boost::numeric::odeint::runge_kutta_dopri5< double > stepper_type;
typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;

using namespace std;

Integrator::Integrator()
{
	printEnergy = false;
	symmetrize = 0;
	integratorType = 0;
	adaptativeError = 1.0e-12;
	defaultTimeStep = 1.0e-3;

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

void Integrator::setOptions(bool printEnergy_in, int symmetrize_in)
{
	//printEnergy = printEnergy_in;
	if (printEnergy)
	{
		printEnergyFile_.open("printEnergyFile.csv");
	}
	symmetrize = symmetrize_in;
}


void Integrator::odeintAdaptativeIntegrator(
	std::vector<double> & xInitial,
	std::vector<double> & vInitial,
	double wholeTimeStep)
{
	FitnessCoulombOdeint function(symmetrize);

	state_type x(xInitial.size() * 2);
	for (int i = 0; i < xInitial.size(); i++)
	{
		x[i] = xInitial[i];
		x[i + 12] = vInitial[i];
	}

	boost::numeric::odeint::integrate_adaptive(
		boost::numeric::odeint::make_controlled< error_stepper_type >
		(adaptativeError, adaptativeError), // error handling
		function,                           // function
		x,                                  // initial state
		0.0e0,                              // initial time
		wholeTimeStep,                      // end_time
		defaultTimeStep);                   // delta time

	for (int i = 0; i < xInitial.size(); i++)
	{
		xInitial[i] = x[i];
		vInitial[i] = x[i + 12];
	}

	// adicionar simetrizacao aqui se for o caso.
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
	
	symm_.symmetrize(symmetrize, xInitial, vInitial);

	/*
	if (symmetrize == 1) // inversion center
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
	else if (symmetrize == 2) // LH2bl (x=0 ; z->-z; y->y)
	{
		xInitial[0] = 0.0e0;
		xInitial[1] = 0.0e0;
		xInitial[2] = 0.0e0;
		xInitial[3] = 0.0e0;
		xInitial[6] = xInitial[4];
		xInitial[7] = xInitial[5];
		xInitial[10] = -xInitial[8];
		xInitial[11] = -xInitial[9];

		vInitial[0] = 0.0e0;
		vInitial[1] = 0.0e0;
		vInitial[2] = 0.0e0;
		vInitial[3] = 0.0e0;

		vInitial[6] = vInitial[4];
		vInitial[7] = vInitial[5];
		vInitial[10] = -vInitial[8];
		vInitial[11] = -vInitial[9];
	}
	else if (symmetrize == 3) // f-langmuirCenter
	{
		xInitial[2] = xInitial[0];
		xInitial[6] = -xInitial[4];
		xInitial[8] = 0.0e0;
		xInitial[10] = 0.0e0;

		vInitial[2] = vInitial[0];
		vInitial[6] = -vInitial[4];
		vInitial[8] = 0.0e0;
		vInitial[10] = 0.0e0;

		xInitial[3] = xInitial[1];
		xInitial[5] = 0.0e0;
		xInitial[7] = 0.0e0;
		xInitial[11] = -xInitial[9];

		vInitial[3] = vInitial[1];
		vInitial[5] = 0.0e0;
		vInitial[7] = 0.0e0;
		vInitial[11] = -vInitial[9];
	}
	else if (symmetrize == 4) // g-langmuirCenter
	{
		xInitial[0] = 0.0e0;
		xInitial[2] = 0.0e0;
		xInitial[6] = -xInitial[4];
		xInitial[10] = xInitial[8];

		vInitial[0] = 0.0e0;
		vInitial[2] = 0.0e0;
		vInitial[6] = -vInitial[4];
		vInitial[10] = vInitial[8];

		xInitial[1] = 0.0e0;
		xInitial[3] = 0.0e0;
		xInitial[5] = 0.0e0;
		xInitial[7] = 0.0e0;

		vInitial[1] = 0.0e0;
		vInitial[3] = 0.0e0;
		vInitial[5] = 0.0e0;
		vInitial[7] = 0.0e0;
	}

	*/
}









