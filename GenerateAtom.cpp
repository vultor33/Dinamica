#include "GenerateAtom.h"

#include <vector>
#include <cmath>

using namespace std;

GenerateAtom::GenerateAtom(int seed)
{
	srand(seed);
	mProton = 1836.15273443449e0;
	pi_ = 3.1415926535897932384626433832795e0;
	
}

GenerateAtom::~GenerateAtom(){}

vector<double> GenerateAtom::addAtom(
	vector<double> & property1,
	vector<double> & property2)
{
	int nParticle1 = property1.size() / 3;
	int nParticle2 = property2.size() / 3;
	vector<double> sum(property1.size() + property2.size());

	for (int i = 0; i < nParticle1; i++)
	{
		sum[i] = property1[i];
		sum[i + nParticle1 + nParticle2] = property1[i + nParticle1];
		sum[i + 2 * (nParticle1 + nParticle2)] = property1[i + 2 * nParticle1];
	}
	for (int i = 0; i < nParticle2; i++)
	{
		sum[i + nParticle1] = property2[i];
		sum[i + nParticle1 + nParticle1 + nParticle2] = property2[i + nParticle2];
		sum[i + nParticle1 + 2 * (nParticle1 + nParticle2)] = property2[i + 2 * nParticle2];
	}

	return sum;
}

vector<double> GenerateAtom::insertVector(vector<double> & vec1, vector<double> &vec2)
{
	vector<double> sum = vec1;
	sum.insert(sum.end(), vec2.begin(), vec2.end());
	return sum;
}

void GenerateAtom::generateTwoRandomAtoms(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	generateInitialPositionAndVelocity(x2, v2, aM2, aC2);
	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);
	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);
}

void GenerateAtom::generateTwoSymmetricAtoms(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	x2 = x1;
	v2 = v1;
	aM2 = aM1;
	aC2 = aC1;
	for (size_t i = 0; i < x1.size(); i++)
	{
		x2[i] *= -1.0e0;
	}
	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);
	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);
}

void GenerateAtom::generateTwoIdenticalSymmetricAtoms(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	x2 = x1;
	v2 = v1;
	aM2 = aM1;
	aC2 = aC1;
	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);
	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);
}


void GenerateAtom::generateTwoAntiSymmetricAtoms(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	x2 = x1;
	v2 = v1;
	aM2 = aM1;
	aC2 = aC1;
	for (size_t i = 0; i < x1.size(); i++)
	{
		x2[i] *= -1.0e0;
		v2[i] *= -1.0e0;
	}
	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);
	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);
}

void GenerateAtom::generateTwoIdenticalAntiSymmetricAtoms(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	x2 = x1;
	v2 = v1;
	aM2 = aM1;
	aC2 = aC1;
	for (size_t i = 0; i < x1.size(); i++)
	{
		v2[i] *= -1.0e0;
	}
	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);
	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);
}


double GenerateAtom::randomNumber(double fMin, double fMax)
{
	double f = ((double)rand() / (double)(RAND_MAX));
	return fMin + f * (fMax - fMin);
}

vector<double> GenerateAtom::unitarySphericalVector()
{
	vector<double> unit(3);
	double fi, teta;
	fi = 2.0e0 * pi_ * randomNumber(0.0e0, 1.0e0);
	teta = acos(2.0e0 * randomNumber(0.0e0, 1.0e0) - 1.0e0);
	unit[0] = sin(teta) * cos(fi);
	unit[1] = sin(teta) * sin(fi);
	unit[2] = cos(teta);
	return unit;
}

void GenerateAtom::translateToCenterOfMass(
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

void GenerateAtom::generateInitialPositionAndVelocity(
	vector<double> &xPositions,
	vector<double> &vVelocities,
	vector<double> &atomsMass,
	vector<double> &atomCharge)
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
	atomCharge.resize(nParticle);
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
	atomCharge[0] = -1.0e0;
	atomCharge[1] = 1.0e0;

	translateToCenterOfMass(xPositions, atomsMass);
}

