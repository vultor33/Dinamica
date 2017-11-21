#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <time.h>

#include "Integrator.h"

using namespace std;

const double mProton = 1836.15273443449e0;

const double pi = 3.1415926535897932384626433832795e0;

void printAtomsVectorDouble(vector<double> & atoms, string testName);

void printCoulombAtoms(vector<double> & atoms, string testName, vector<double> &atomsCharge);

void translateToCenterOfMass(vector<double> &x, vector<double> &atomsMass);

vector<double> unitarySphericalVector();

double randomNumber(double fMin, double fMax);

void generateInitialPositionAndVelocity(
	vector<double> &xPositions, 
	vector<double> &vVelocities,
	vector<double> &atomsMass);



int main()
{
//	int seed = 3;
	srand(time(NULL));

	// SCF = wait until electrons find an orbit

	/*
	int natm = 2;
	vector<double> x(3 * natm);
	vector<double> v(3 * natm);
	vector<double> atomsMass(3 * natm);
	vector<double> atomsCharge(natm);
	atomsCharge[0] = 1.0e0;
	atomsCharge[1] = -1.0e0;
	for (size_t i = 0; i < x.size(); i++)
	{
		x[i] = 0.0e0;
		v[i] = 0.0e0;
		atomsMass[i] = 1.0e0;
	}
	//1
	x[0] = 0.0e0;
	atomsMass[0] = mProton;
	atomsMass[0 + natm] = mProton;
	atomsMass[0 + 2 * natm] = mProton;

	//2
	x[1] = 1.0e0;
	v[1 + natm] = 1.0e0;

	//3
//	x[5] = 1;
//	x[8] = 0.4;
*/

	vector<double> x1;
	vector<double> v1;
	vector<double> atomsMass1;
	vector<double> atomsCharge1(2);
	atomsCharge1[0] = -1.0e0;
	atomsCharge1[1] = 1.0e0;
	generateInitialPositionAndVelocity(x1, v1, atomsMass1);
	vector<double> x2;
	vector<double> v2;
	vector<double> atomsMass2;
	vector<double> atomsCharge2(2);
	atomsCharge2[0] = -1.0e0;
	atomsCharge2[1] = 1.0e0;
	generateInitialPositionAndVelocity(x2, v2, atomsMass2);

	// translade cm1 -10x e 1y.
	x1[0] -= 5.0e0;
	x1[1] -= 5.0e0;
	x1[2] += 2.0e0;
	x1[3] += 2.0e0;
	v1[0] += 0.001e0;
	v1[1] += 0.001e0;
	v2[0] -= 0.001e0;
	v2[1] -= 0.001e0;

	printCoulombAtoms(x1, "printInitial.xyz", atomsCharge1);
	printCoulombAtoms(x2, "printInitial.xyz", atomsCharge2);

	vector<double> x(12);
	vector<double> v(12);
	vector<double> atomsMass(12);
	vector<double> atomsCharge(4);
	x[0] = x1[0];
	x[1] = x1[1];
	x[2] = x2[0];
	x[3] = x2[1];
	x[4] = x1[2];
	x[5] = x1[3];
	x[6] = x2[2];
	x[7] = x2[3];
	x[8] = x1[4];
	x[9] = x1[5];
	x[10] = x2[4];
	x[11] = x2[5];
	v[0] = v1[0];
	v[1] = v1[1];
	v[2] = v2[0];
	v[3] = v2[1];
	v[4] = v1[2];
	v[5] = v1[3];
	v[6] = v2[2];
	v[7] = v2[3];
	v[8] = v1[4];
	v[9] = v1[5];
	v[10] = v2[4];
	v[11] = v2[5];
	atomsMass[0] = atomsMass1[0];
	atomsMass[1] = atomsMass1[1];
	atomsMass[2] = atomsMass2[0];
	atomsMass[3] = atomsMass2[1];
	atomsMass[4] = atomsMass1[2];
	atomsMass[5] = atomsMass1[3];
	atomsMass[6] = atomsMass2[2];
	atomsMass[7] = atomsMass2[3];
	atomsMass[8] = atomsMass1[4];
	atomsMass[9] = atomsMass1[5];
	atomsMass[10] = atomsMass2[4];
	atomsMass[11] = atomsMass2[5];
	atomsCharge[0] = atomsCharge1[0];
	atomsCharge[1] = atomsCharge1[1];
	atomsCharge[2] = atomsCharge2[0];
	atomsCharge[3] = atomsCharge2[1];

	printCoulombAtoms(x, "printInitial.xyz", atomsCharge);

	translateToCenterOfMass(x, atomsMass);

	printCoulombAtoms(x, "printInitial.xyz", atomsCharge);

	Integrator rk_;
	rk_.setAdditionalParams(atomsMass, atomsCharge);
	rk_.setOptions(true);
	remove("simulacao.xyz");
	printCoulombAtoms(x, "simulacao.xyz", atomsCharge);
	int iterationsLoop = 10;
	int printLoop = 15000;
	int maxIterations = printLoop * iterationsLoop;
	for (int i = 0; i < printLoop;i++)
	{
		for (int j = 0; j < iterationsLoop; j++)
		{
			rk_.rungeKuttaSimetrico(x, v, 0.01);
		}
		cout << 100 * i / printLoop << " %" << endl; ;
		printCoulombAtoms(x, "simulacao.xyz", atomsCharge);
		for (int i = 0; i < 12; i++)
		{
			if ((x[i] > 7.0e0) || (x[i] < -7.0e0))
				v[i] *= -1.0e0;
		}
	}

	return 0;


}

double randomNumber(double fMin, double fMax)
{
	double f = ((double)rand() / (double)(RAND_MAX));
	return fMin + f * (fMax - fMin);
}

vector<double> unitarySphericalVector()
{
	vector<double> unit(3);
	double fi, teta;
	fi = 2.0e0 * pi * randomNumber(0.0e0, 1.0e0);
	teta = acos(2.0e0 * randomNumber(0.0e0, 1.0e0) - 1.0e0);
	unit[0] = sin(teta) * cos(fi);
	unit[1] = sin(teta) * sin(fi);
	unit[2] = cos(teta);
	return unit;
}

void translateToCenterOfMass(
	vector<double> &x, 
	vector<double> &atomsMass)
{
	double cmx = 0.0e0;
	double cmy = 0.0e0;
	double cmz = 0.0e0;
	int natm = x.size() / 3;
	double totalMass = 0.0e0;
	for (int i = 0; i < natm; i++)
	{
		totalMass += atomsMass[i];
		cmx += x[i] * atomsMass[i];
		cmy += x[i + natm] * atomsMass[i + natm];
		cmz += x[i + 2 * natm] * atomsMass[i + 2 * natm];
	}
	cmx /= totalMass;
	cmy /= totalMass;
	cmz /= totalMass;
	for (int i = 0; i < natm; i++)
	{
		x[i] -= cmx;
		x[i + natm] -= cmy;
		x[i + 2 * natm] -= cmz;
	}
}

void generateInitialPositionAndVelocity(
	vector<double> &xPositions, 
	vector<double> &vVelocities,
	vector<double> &atomsMass)
{
	vector<double> ePos = unitarySphericalVector();
	vector<double> eVel = unitarySphericalVector();
	double prodIntVelPos = eVel[0] * ePos[0] + eVel[1] * ePos[1] + eVel[2] * ePos[2];
	// project on perpendicular plane
	eVel[0] = eVel[0] - prodIntVelPos * ePos[0];
	eVel[1] = eVel[1] - prodIntVelPos * ePos[1];
	eVel[2] = eVel[2] - prodIntVelPos * ePos[2];

	double norm = sqrt(eVel[0] * eVel[0] + eVel[1] * eVel[1] + eVel[2] * eVel[2]);
	eVel[0] /= norm;
	eVel[1] /= norm;
	eVel[2] /= norm;

	double perpend = eVel[0] * ePos[0] + eVel[1] * ePos[1] + eVel[2] * ePos[2];

	vector<double> protonVel = eVel;
	protonVel[0] /= -mProton;
	protonVel[1] /= -mProton;
	protonVel[2] /= -mProton;

	int nParticle = 2;
	xPositions.resize(nParticle * 3);
	vVelocities.resize(nParticle * 3);
	atomsMass.resize(nParticle * 3);
	xPositions[0] = ePos[0];
	xPositions[nParticle] = ePos[1];
	xPositions[2 * nParticle] = ePos[2];
	xPositions[1] = 0.0e0;
	xPositions[1 + nParticle] = 0.0e0;
	xPositions[1 + 2 * nParticle] = 0.0e0;
	vVelocities[0] = eVel[0];
	vVelocities[nParticle] = eVel[1];
	vVelocities[2 * nParticle] = eVel[2];
	vVelocities[1] = protonVel[0];
	vVelocities[1 + nParticle] = protonVel[1];
	vVelocities[1 + 2 * nParticle] = protonVel[2];
	atomsMass[0] = 1.0e0;
	atomsMass[nParticle] = 1.0e0;
	atomsMass[2 * nParticle] = 1.0e0;
	atomsMass[1] = mProton;
	atomsMass[1 + nParticle] = mProton;
	atomsMass[1 + 2 * nParticle] = mProton;

	translateToCenterOfMass(xPositions, atomsMass);
}


void printAtomsVectorDouble(vector<double> & atoms, string testName)
{
	int natm = atoms.size() / 3;
	ofstream teste_;
	teste_.open(testName.c_str(), std::ofstream::out | std::ofstream::app);
	teste_ << natm << endl << "t" << endl;
	for (int i = 0; i < natm; i++)
	{
		teste_ << "H "
			<< atoms[i] << "  "
			<< atoms[i + natm] << "  "
			<< atoms[i + 2 * natm] << endl;
	}
	teste_.close();
}

void printCoulombAtoms(vector<double> & atoms, string testName, vector<double> &atomsCharge)
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



/*
DIDNT WORKED
double normV1 = sqrt(v[0] * v[0] + v[4] * v[4] + v[8] * v[8]);
if (normV1 > 10.0e0)
{
v[0] *= 0.99;
v[4] *= 0.99;
v[8] *= 0.99;
}
double normV2 = sqrt(v[2] * v[2] + v[6] * v[6] + v[10] * v[10]);
if (normV2 > 10.0e0)
{
v[2] *= 0.99;
v[6] *= 0.99;
v[10] *= 0.99;
}

*/
