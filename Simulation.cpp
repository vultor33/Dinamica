#include "Simulation.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "GenerateAtom.h"
#include "Integrator.h"
#include "DynamicsStructs.h"

using namespace std;

Simulation::Simulation(DymOptions & dymOptions_in)
{
	dymOptions_ = dymOptions_in;

	stopSimulation = false;

}

Simulation::~Simulation(){}

void Simulation::startSimulation(
	int seed, 
	double tempKelvin, 
	double impactFactorAu, 
	int simulationType
)
{
	GenerateAtom genAtom_(dymOptions_.seed);
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
		dymOptions_.symmetrize = true;
		break;

	case 6:
		genAtom_.generateBohrMolecule(x, v, atomsMass, atomsCharge);
		dymOptions_.symmetrize = true;
		break;

	default:
		cout << "simulation type not found" << endl;
		exit(1);
		break;
	}
	x[0] -= dymOptions_.initialDistance;
	x[1] -= dymOptions_.initialDistance;
	x[4] += dymOptions_.impactParameter;
	x[5] += dymOptions_.impactParameter;
	v[0] += dymOptions_.initialSpeed;
	v[1] += dymOptions_.initialSpeed;
	v[2] -= dymOptions_.initialSpeed;
	v[3] -= dymOptions_.initialSpeed;
	genAtom_.translateToCenterOfMass(x, atomsMass);

	Integrator rk_;
	rk_.setAdditionalParams(atomsMass, atomsCharge);
	rk_.setOptions(dymOptions_.printEnergy, dymOptions_.symmetrize);

	if (dymOptions_.printMovie)
	{
		remove(dymOptions_.outName.c_str());
		printCoulombAtoms(x, dymOptions_.outName, atomsCharge);
	}
	if (dymOptions_.printPosVel)
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
	if (dymOptions_.printPosVel)
		printPositionsAndVelocities(x, v, posVel_);
	for (int i = 0; i < dymOptions_.printLoop; i++)
	{
		for (int j = 0; j < dymOptions_.iterationLoop; j++)
		{
			rk_.rungeKuttaSimetrico(x, v, dymOptions_.timeStep);
		}
		if (dymOptions_.checkStopSimulationConditions)
			checkStopSimulation(x);
		if (stopSimulation)
			break;
		if (dymOptions_.printPosVel)
			printPositionsAndVelocities(x, v, posVel_);

		if (dymOptions_.printMovie)
		{
			cout << 100 * i / dymOptions_.printLoop << " %" << endl; ;
			printCoulombAtoms(x, dymOptions_.outName, atomsCharge);
		}
	}

}

void Simulation::additionalOptions(string flag, bool option)
{
	if (flag == "printMovie")
		dymOptions_.printMovie = option;
	else if (flag == "printEnergy")
		dymOptions_.printEnergy = option;
	else if (flag == "printPosVel")
		dymOptions_.printPosVel = option;
	else
	{
		cout << "flag not found" << endl;
		exit(1);
	}
}


void Simulation::checkStopSimulation(vector<double> & x)
{
	for (size_t i = 0; i < x.size(); i++)
	{
		if (abs(x[i]) > dymOptions_.maxStopSimulationDistance)
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

	/* MOMENTO ANGULAR
	LxE1 =(I2*R2-Q2*J2)
	LyE1 =(Q2*B2-A2*R2)
	LzE1 =A2*J2-I2*B2
	LxE2 =M2*V2-U2*N2
	LyE2 =U2*F2-E2*V2
	LzE2 =E2*N2-M2*F2
	LxP1 =(K2*T2-S2*L2)*1836.15273443449
	LyP1 =(S2*D2-C2*T2)*1836.15273443449
	LzP1 =(C2*L2-K2*D2)*1836.15273443449
	LxP2 =(O2*X2-W2*P2)*1836.15273443449
	LyP2 =(W2*H2-X2*G2)*1836.15273443449
	LzP2 =(G2*P2-O2*H2)*1836.15273443449

	LxTOTAL =AJ2+AG2+AD2+AA2
	LyTOTAL =AK2+AH2+AE2+AB2
	LzTOTAL =AL2+AI2+AF2+AC2

	L =RAIZ(AN2*AN2+AO2*AO2+AP2*AP2)

	*/



}
