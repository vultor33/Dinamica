#include "Simulation.h"

#include <vector>
#include <iostream>

#include "GenerateAtom.h"
#include "Integrator.h"

using namespace std;

Simulation::Simulation()
{
	timeStep = 0.01;
	iterationLoop = 100;
	printLoop = 5000;
	initialDistance = 5.0e0;
	impactParameter = 1.0e0;
	initialSpeed = 0.0005;
	checkStopSimulationConditions = false;
	stopSimulation = false;
	printEnergy = true;
	simmetrize = false;
	outputName = "simulacao.xyz";

	// CENTRO DE MASSA MEXENDO - TEM MAIS COISA AQUI.



}

Simulation::~Simulation(){}

void Simulation::startSimulation(int seed)
{
	GenerateAtom genAtom_(seed);
	vector<double> x, v, atomsMass, atomsCharge;
	genAtom_.generateTwoAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
	x[0] -= initialDistance;
	x[1] -= initialDistance;
	x[4] += impactParameter;
	x[5] += impactParameter;
	v[0] += initialSpeed;
	v[1] += initialSpeed;
	v[2] -= initialSpeed;
	v[3] -= initialSpeed;
	genAtom_.translateToCenterOfMass(x, atomsMass);

	Integrator rk_;
	rk_.setAdditionalParams(atomsMass, atomsCharge);
	rk_.setOptions(printEnergy, simmetrize);

	remove(outputName.c_str());
	printCoulombAtoms(x, outputName, atomsCharge);
	int maxIterations = printLoop * iterationLoop;
	for (int i = 0; i < printLoop; i++)
	{
		for (int j = 0; j < iterationLoop; j++)
		{
			rk_.rungeKuttaSimetrico(x, v, timeStep);
		}
		cout << 100 * i / printLoop << " %" << endl; ;
		printCoulombAtoms(x, outputName, atomsCharge);
		if (checkStopSimulationConditions)
			checkStopSimulation(x);
		if (stopSimulation)
			break;
	}

}


void Simulation::checkStopSimulation(vector<double> & x)
{
	for (size_t i = 0; i < x.size(); i++)
	{
		if (abs(x[i]) > 20.0e0)
			stopSimulation = true;
	}
}

void Simulation::printCoulombAtoms(vector<double> & atoms, string testName, vector<double> &atomsCharge)
{
	int natm = atoms.size() / 3;
	ofstream teste_;
	teste_.open(testName.c_str(), std::ofstream::out | std::ofstream::app);
	teste_ << natm << endl << "t" << endl;
	for (int i = 0; i < natm; i++)
	{
		if (atomsCharge[i] == 1.0e0)
			teste_ << "Au ";
		else
			teste_ << "H ";

		teste_ << atoms[i] << "  "
			<< atoms[i + natm] << "  "
			<< atoms[i + 2 * natm] << endl;
	}
	teste_.close();
}
