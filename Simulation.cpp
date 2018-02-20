#include "Simulation.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>

#include "GenerateAtom.h"
#include "Integrator.h"

using namespace std;

Simulation::Simulation()
{
	timeStep = 0.01e0;
	// 5000000 it * print
	iterationLoop = 1000;
	printLoop = 300;
	//initialDistance = 5.0e0;
	initialDistance = 0.0e0;
	impactParameter = 1.0e0;
	initialSpeed = 0.001;
	checkStopSimulationConditions = true;
	stopSimulation = false;
	printEnergy = false;
	simmetrize = false;
	printMovie = false;
	printPosVel = true;
	maxStopSimulationDistance = 100.0e0;
	temperatureUnit = 315774.64e0;
	mProton = 1836.15273443449e0;
	outputName = "simulacao.xyz";
	// CENTRO DE MASSA MEXENDO UM POUCO (10^-13)
}

Simulation::~Simulation(){}

void Simulation::startSimulation(
	int seed, 
	double tempKelvin, 
	double impactFactorAu, 
	int simulationType
)
{
	initialVelocityKinecticTheory(tempKelvin);
	impactParameter = impactFactorAu;
	GenerateAtom genAtom_(seed);
	vector<double> x, v, atomsMass, atomsCharge;
	switch (simulationType)
	{
	case 0:
		genAtom_.generateTwoRandomAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 1:
		genAtom_.generateTwoAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 2:
		genAtom_.generateTwoIdenticalAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 3:
		genAtom_.generateTwoSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 4:
		genAtom_.generateTwoIdenticalSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 5:
		genAtom_.generateTwoAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
		simmetrize = true;
		break;

	case 6:
		genAtom_.generateBohrMolecule(x, v, atomsMass, atomsCharge);
		simmetrize = true;
		break;
		

	default:
		cout << "simulation type not found" << endl;
		exit(1);
		break;
	}
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

	if (printMovie)
	{
		remove(outputName.c_str());
		printCoulombAtoms(x, outputName, atomsCharge);
	}
	if (printPosVel)
	{
		stringstream buildName;
		buildName << seed << "-"
			<< tempKelvin << "-"
			<< impactFactorAu << "-"
			<< simulationType;

		string posVelName = "simulation-" + buildName.str() + ".csv";
		posVel_.open(posVelName);
		posVel_ << "xE1 ; vxE1 ; xP1 ; vxP1 ; xE2 ; vxE2; xP2 ; vxP2 ;"
			<< " yE1 ; vyE1 ; yP1 ; vyP1 ; yE2 ; vyE2; yP2 ; vyP2 ;"
			<< " zE1 ; vzE1 ; zP1 ; vzP1 ; zE2 ; vzE2; zP2 ; vzP2 "
			<< endl;
	}
	if (printPosVel)
		printPositionsAndVelocities(x, v, posVel_);
	for (int i = 0; i < printLoop; i++)
	{
		for (int j = 0; j < iterationLoop; j++)
		{
			rk_.rungeKuttaSimetrico(x, v, timeStep);
		}
		if (checkStopSimulationConditions)
			checkStopSimulation(x);
		if (stopSimulation)
			break;
		if (printPosVel)
			printPositionsAndVelocities(x, v, posVel_);

		if (printMovie)
		{
			cout << 100 * i / printLoop << " %" << endl; ;
			printCoulombAtoms(x, outputName, atomsCharge);
		}
	}

}

void Simulation::additionalOptions(string flag, bool option)
{
	if (flag == "printMovie")
		printMovie = option;
	else if (flag == "printEnergy")
		printEnergy = option;
	else if (flag == "printPosVel")
		printPosVel = option;
	else
	{
		cout << "flag not found" << endl;
		exit(1);
	}
}

void Simulation::initialVelocityKinecticTheory(double TempKelvin)
{
	// v = sqrt( 3 k t / m)
	// i am suggested this half factor
	double TempAUnits = TempKelvin / temperatureUnit;
	initialSpeed = sqrt(3.0e0 * TempAUnits / (1.0e0 + mProton)) / 2.0e0;
}
void Simulation::checkStopSimulation(vector<double> & x)
{
	for (size_t i = 0; i < x.size(); i++)
	{
		if (abs(x[i]) > maxStopSimulationDistance)
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

void Simulation::printPositionsAndVelocities(
	std::vector<double> &x,
	std::vector<double> &v,
	ofstream &posVelFile_)
{
	for (size_t i = 0; i < x.size(); i++)
	{
		posVelFile_ << x[i] << " ; " << v[i] << " ; ";
	}
	posVelFile_ << endl;
}
