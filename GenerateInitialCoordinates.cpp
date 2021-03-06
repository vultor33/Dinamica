#include "GenerateInitialCoordinates.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "DynamicsStructs.h"
#include "AuxMath.h"
#include "TrajectorySymmetrizer.h"

using namespace std;

GenerateInitialCoordinates::GenerateInitialCoordinates()
{
	mProton = 1836.15273443449e0;
	pi_ = 3.1415926535897932384626433832795e0;
}

GenerateInitialCoordinates::~GenerateInitialCoordinates(){}

bool GenerateInitialCoordinates::generateInitial(
	DymOptions &dymOptions_,
	std::vector<double> &x,
	std::vector<double> &v,
	std::vector<double> &atomsMass,
	std::vector<double> &atomsCharge)
{
	bool sucess = true;

	switch (dymOptions_.simulationType)
	{
	case 0:
		generateTwoRandomAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 1:
		generateTwoAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 2:
		generateTwoIdenticalAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 3:
		generateTwoSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 4:
		generateTwoIdenticalSymmetricAtoms(x, v, atomsMass, atomsCharge);
		break;

	case 5:
		generateTwoAntiSymmetricAtoms(x, v, atomsMass, atomsCharge);
		dymOptions_.symmetrize = 1;
		break;

	case 6:
		generateBohrMolecule(x, v, atomsMass, atomsCharge, dymOptions_.energy, dymOptions_.angle);
		break;

	case 7:
		generateBohrEllipseMolecule(
			x, 
			v, 
			atomsMass, 
			atomsCharge, 
			dymOptions_.energy, 
			dymOptions_.rElec,
			dymOptions_.rProton);
		break;

	case 8:
		generateBohrEllipseMoleculeAngle(
			x,
			v,
			atomsMass,
			atomsCharge,
			dymOptions_.energy,
			dymOptions_.rElec,
			dymOptions_.rProton,
			dymOptions_.angle);
		break;
	case 9:
		generateBohrEllipseMoleculeAngle(
			x,
			v,
			atomsMass,
			atomsCharge,
			dymOptions_.energy,
			dymOptions_.rElec,
			dymOptions_.rProton,
			180.0e0);
		break;

	case 10:
		sucess = generateElectronsAtCenter(
			x,
			v,
			atomsMass,
			atomsCharge,
			dymOptions_);
		if (!sucess)
			return false;
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
	translateToCenterOfMass(x, atomsMass);
	return sucess;

}


vector<double> GenerateInitialCoordinates::addAtom(
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

vector<double> GenerateInitialCoordinates::insertVector(vector<double> & vec1, vector<double> &vec2)
{
	vector<double> sum = vec1;
	sum.insert(sum.end(), vec2.begin(), vec2.end());
	return sum;
}

void GenerateInitialCoordinates::generateTwoRandomAtoms(
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

void GenerateInitialCoordinates::generateTwoSymmetricAtoms(
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

void GenerateInitialCoordinates::generateTwoIdenticalSymmetricAtoms(
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


void GenerateInitialCoordinates::generateTwoAntiSymmetricAtoms(
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

void GenerateInitialCoordinates::generateTwoIdenticalAntiSymmetricAtoms(
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


void GenerateInitialCoordinates::generateBohrMolecule(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge,
	double energy,
	double angle)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	// x[0, 2 e 4] eletron
	// x[1, 3 e 5] proton
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	for (size_t i = 0; i < x1.size(); i++)
	{
		x1[i] = 0.0e0;
		v1[i] = 0.0e0;
	}

	double rProton, rElec, vInit;
	calcBohrParams(energy, rProton, rElec, vInit);

	x1[2] = rElec;
	x1[5] = rProton;

	AuxMath auxMath_;
	double vx, vz;
	if (abs(angle - 90) < 0.00001)
	{
		vx = 0.0e0;
		vz = vInit;
	}
	else
	{
		double tanAngle = tan(angle * auxMath_._pi / (180.0e0));
		vx = sqrt(vInit*vInit / (1 + tanAngle * tanAngle));
		vz = vx * tanAngle;
	}

	v1[0] = vx;
	v1[4] = vz;

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

void GenerateInitialCoordinates::calcBohrParams(
	double energy,
	double &rProtonInit,
	double &rElecInit,
	double &vInit)
{
	rProtonInit = (3.0e0 - 9.0e0*sqrt(3.0e0)) / (12.0e0*sqrt(3.0e0)*energy);

	rElecInit = rProtonInit * sqrt(3.0e0);

	vInit = sqrt((9.0e0 - sqrt(3.0e0)) / (12.0e0 * rProtonInit));

}


void GenerateInitialCoordinates::generateBohrEllipseMoleculeAngle(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge,
	double energy,
	double rElec,
	double rProton,
	double angle)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	// x[0, 2 e 4] eletron
	// x[1, 3 e 5] proton
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	for (size_t i = 0; i < x1.size(); i++)
	{
		x1[i] = 0.0e0;
		v1[i] = 0.0e0;
	}

	x1[4] = rElec;
	x1[5] = rProton;

	double vPotential = (1.0e0 / (2.0e0*rProton)) 
		+ (1.0e0 / (2.0e0*rElec)) 
		- (2.0e0 / (rElec - rProton)) 
		- (2.0e0 / (rElec + rProton));

	if ((energy - vPotential) < 0.0e0)
		v1[2] = 10000.0e0;
	else
		v1[2] = sqrt(energy - vPotential);

	x2 = x1;
	v2 = v1;
	aM2 = aM1;
	aC2 = aC1;
	for (size_t i = 0; i < x1.size(); i++)
	{
		x2[i] *= -1.0e0;
		v2[i] *= -1.0e0;
	}

	double vMod = v2[2];
	AuxMath auxMath_;
	double vx, vy;
	if ((abs(angle - 90) < 0.00001))
	{
		vy = 0.0e0;
		vx = vMod;
	}
	else if (abs(angle - 90) < 0.00001)
	{
		vy = 0.0e0;
		vx = -vMod;
	}
	else
	{
		double tanAngle = tan(angle * auxMath_._pi / (180.0e0));
		vy = sqrt(vMod*vMod / (1 + tanAngle * tanAngle));
		vx = vy * tanAngle;
	}
	v2[2] = vy;
	v2[0] = vx;


	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);
	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);
}

void GenerateInitialCoordinates::generateBohrEllipseMolecule(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge,
	double energy,
	double rElec,
	double rProton)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	// x[0, 2 e 4] eletron
	// x[1, 3 e 5] proton
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	for (size_t i = 0; i < x1.size(); i++)
	{
		x1[i] = 0.0e0;
		v1[i] = 0.0e0;
	}

	x1[4] = rElec;
	x1[5] = rProton;

	double vPotential = (1.0e0 / (2.0e0*rProton))
		+ (1.0e0 / (2.0e0*rElec))
		- (2.0e0 / (rElec - rProton))
		- (2.0e0 / (rElec + rProton));

	if ((energy - vPotential) < 0.0e0)
		v1[2] = 10000.0e0;
	else
		v1[2] = sqrt(energy - vPotential);

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



double GenerateInitialCoordinates::randomNumber(double fMin, double fMax)
{
	double f = ((double)rand() / (double)(RAND_MAX));
	return fMin + f * (fMax - fMin);
}

vector<double> GenerateInitialCoordinates::unitarySphericalVector()
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

void GenerateInitialCoordinates::translateToCenterOfMass(
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

void GenerateInitialCoordinates::velocityCmCorrection(std::vector<double> &v, std::vector<double> &atomsMass)
{
	double vmx = 0.0e0;
	double vmy = 0.0e0;
	double vmz = 0.0e0;
	int natm = v.size() / 3;
	double mProton = atomsMass[1];
	for (int i = 0; i < natm; i++)
	{
		vmx += v[i] * atomsMass[i];
		vmy += v[i + natm] * atomsMass[i + natm];
		vmz += v[i + 2 * natm] * atomsMass[i + 2 * natm];
	}
	double vPCorrX = vmx / (2.0e0 + 2.0e0 * mProton);
	double vPCorrY = vmy / (2.0e0 + 2.0e0 * mProton);
	double vPCorrZ = vmz / (2.0e0 + 2.0e0 * mProton);
	v[1] -= vPCorrX;
	v[3] -= vPCorrX;
	v[5] -= vPCorrY;
	v[7] -= vPCorrY;
	v[9] -= vPCorrZ;
	v[11] -= vPCorrZ;
	vmx = vmy = vmz = 0.0e0;
	for (int i = 0; i < natm; i++)
	{
		vmx += v[i] * atomsMass[i];
		vmy += v[i + natm] * atomsMass[i + natm];
		vmz += v[i + 2 * natm] * atomsMass[i + 2 * natm];
	}
	vmx = 1;
}

// generate 1 eletron and 1 proton. 
// x[0, 2 e 4] eletron
// x[1, 3 e 5] proton
void GenerateInitialCoordinates::generateInitialPositionAndVelocity(
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


























// DEFINICAO DO NOVO CONJUNTO DE TRAJETORIAS
bool GenerateInitialCoordinates::generateElectronsAtCenter(
	std::vector<double> &xPositions,
	std::vector<double> &vVelocities,
	std::vector<double> &atomsMass,
	std::vector<double> &atomCharge,
	DymOptions &dymOptions_)
{
	vector<double> x1, x2, v1, v2, aM1, aM2, aC1, aC2;
	generateInitialPositionAndVelocity(x1, v1, aM1, aC1);
	for (size_t i = 0; i < x1.size(); i++)
	{
		x1[i] = 0.0e0;
		v1[i] = 0.0e0;
	}

	double energy = dymOptions_.energy;
	double rElec = dymOptions_.rElec;
	double rProton = dymOptions_.rProton;
	double angle = dymOptions_.angle;

	double vPotential;
	if (dymOptions_.initialPositionType == 0)
	{
		vPotential = (1.0e0 / (2.0e0*rProton))
			+ (1.0e0 / (2.0e0*rElec))
			- (4.0e0 / sqrt(rElec*rElec + rProton * rProton));
	}
	else if (dymOptions_.initialPositionType == 1)
	{
		vPotential = (1.0e0 / (2.0e0*rProton))
			+ (1.0e0 / (2.0e0*rElec))
			- (2.0e0 / (rElec - rProton))
			- (2.0e0 / (rElec + rProton));
	}

	double vInitial;
	if ((energy - vPotential) < 0.0e0)
	{
		return false;
	}

	if (dymOptions_.initialPositionType == 0)
	{
		x1[2] = rElec; //ye
		x1[5] = rProton; //zp
	}
	else if (dymOptions_.initialPositionType == 1)
	{
		x1[4] = rElec; //ze
		x1[5] = rProton; //zp
	}

	/*
	vInitial = sqrt(energy - vPotential);

	AuxMath auxMath_;
	double auxV1, auxV2;
	if (abs(angle - 90) < 0.00001)
	{
		auxV1 = 0.0e0;
		auxV2 = vInitial;
	}
	else
	{
		double tanAngle = tan(angle * auxMath_._pi / (180.0e0));
		auxV1 = sqrt(vInitial*vInitial / (1 + tanAngle * tanAngle));
		auxV2 = auxV1 * tanAngle;
	}

	if (dymOptions_.initialPositionType == 0)
	{
		v1[0] = auxV1; //vxe
		v1[4] = auxV2; //vze
	}
	else if (dymOptions_.initialPositionType == 1)
	{
		v1[2] = auxV1; //vye
		v1[0] = auxV2; //vxe
	}
	*/

	//generating initial velocities with a correction at the center of mass
	initialVelocities(energy - vPotential, v1, dymOptions_);

	x2 = x1;
	v2 = v1;
	aM2 = aM1;
	aC2 = aC1;

	xPositions = addAtom(x1, x2);
	vVelocities = addAtom(v1, v2);

	if (dymOptions_.initialPositionType == 0)
	{
		xPositions[6] = -xPositions[4];
		xPositions[11] = -xPositions[9];
		vVelocities[6] = vVelocities[2];
		vVelocities[11] = vVelocities[9];
	}
	else if (dymOptions_.initialPositionType == 1)
	{
		xPositions[10] = -xPositions[8];
		xPositions[11] = -xPositions[9];
		vVelocities[10] = vVelocities[8];
		vVelocities[11] = vVelocities[9];
	}

	TrajectorySymmetrizer symm_;
	symm_.symmetrize(dymOptions_.symmetrize, xPositions, vVelocities);

	atomsMass = addAtom(aM1, aM2);
	atomCharge = insertVector(aC1, aC2);

	return true;
}


void GenerateInitialCoordinates::initialVelocities(
	double deltaEnergy,
	vector<double> &v,
	DymOptions &dymOptions_)
{
	// duas condicoes - centro e canto.


	double vInitial;
	AuxMath auxMath_;
	double auxV1, auxV2;
	double tanTeta, tanSquare, mOneTanSqrt;

	if (abs(dymOptions_.angle - 90) > 0.00001)
	{
		tanTeta = tan(auxMath_._pi * dymOptions_.angle / 180.0e0);
		tanSquare = tanTeta * tanTeta;
		mOneTanSqrt = mProton * (1.0e0 + tanSquare);
	}

	if (dymOptions_.initialPositionType == 0)
	{
		if (dymOptions_.symmetrize == 1)
		{
			vInitial = sqrt(deltaEnergy);
			if (abs(dymOptions_.angle - 90) < 0.00001)
			{
				auxV1 = 0.0e0;
				auxV2 = vInitial;
			}
			else
			{
				auxV1 = vInitial * sqrt(1.0e0 / (1 + tanSquare));
				auxV2 = auxV1 * tanTeta;
			}
		}
		if (dymOptions_.symmetrize == 2)
		{
			vInitial = sqrt((mOneTanSqrt / (mOneTanSqrt + 1.0e0)) * deltaEnergy);
			auxV1 = vInitial * sqrt(1.0e0 / (1 + tanSquare));
			auxV2 = auxV1 * tanTeta;
			v[1] = -auxV1 / mProton; //vpe
		}
		else if (dymOptions_.symmetrize == 3)
		{
			if (abs(dymOptions_.angle - 90) < 0.00001)
			{
				vInitial = sqrt((mProton / (mProton + 1.0e0)) * deltaEnergy);
				auxV1 = 0.0e0;
				auxV2 = vInitial;
				v[5] = -auxV2 / mProton; //vpe
			}
			else
			{
				vInitial = sqrt((mOneTanSqrt / (mOneTanSqrt + tanSquare)) * deltaEnergy);
				auxV1 = vInitial * sqrt(1.0e0 / (1 + tanSquare));
				auxV2 = auxV1 * tanTeta;
				v[5] = -auxV2 / mProton; //vpe
			}
		}
		else if (dymOptions_.symmetrize == 4)
		{
			vInitial = sqrt((mProton / (mProton + 1.0e0)) * deltaEnergy);
			auxV1 = vInitial * sqrt(1.0e0 / (1 + tanSquare));
			auxV2 = auxV1 * tanTeta;
			v[1] = -auxV1 / mProton;
			v[5] = -auxV2 / mProton; //vpe
		}
		v[0] = auxV1; //vxe
		v[4] = auxV2; //vze
	}
	else if (dymOptions_.initialPositionType == 1)
	{
		if (dymOptions_.symmetrize == 1)
		{
			vInitial = sqrt(deltaEnergy);
			if (abs(dymOptions_.angle - 90) < 0.00001)
			{
				auxV1 = 0.0e0;
				auxV2 = vInitial;
			}
			else
			{
				auxV1 = vInitial * sqrt(1.0e0 / (1 + tanSquare));
				auxV2 = auxV1 * tanTeta;
			}
		}
		else if (dymOptions_.symmetrize == 2)
		{
			vInitial = sqrt((mOneTanSqrt / (mOneTanSqrt + 1.0e0)) * deltaEnergy);
			auxV1 = vInitial * sqrt(1.0e0 / (1 + tanSquare));
			auxV2 = auxV1 * tanTeta;
			v[1] = -auxV1 / mProton;
		}
		v[0] = auxV1; //vxe
		v[2] = auxV2; //vye
	}

}

