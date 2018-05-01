#include "Fitness.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include "AuxMath.h"

using namespace std;

Fitness::Fitness()
{
	atomsChargeFit.resize(4);
	mProton = 1836.15273443449e0;
	atomsChargeFit[0] = -1.0e0;
	atomsChargeFit[1] = 1.0e0;
	atomsChargeFit[2] = -1.0e0;
	atomsChargeFit[3] = 1.0e0;
}

Fitness::~Fitness(){}

double Fitness::fit(vector<double> &point, int type)
{
	switch (type)
	{
	case 0:
		return lennardJones(point);
		break;
		
	default:
		cout << "FITNESS FUNCTION NOT FOUND" << endl;
		exit(3);
	}
}

double Fitness::lennardJones(vector<double> &x)
{
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	double r, r2, r4, r6, r12;
	double vlj = 0.0e0;
	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
				);
			r2 = r * r;
			r4 = r2 * r2;
			r6 = r4 * r2;
			r12 = r6 * r6;
			vlj += 4.0e0 * (-1 / r6 + 1 / r12);
		}
	}
	if (isnan(vlj))
		return 1.0e99;
	return vlj;
}


bool Fitness::lennardJonesGradient(vector<double> &x, vector<double> &gradient)
{
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	for (size_t i = 0; i < gradient.size(); i++)
		gradient[i] = 0.0e0;

	double r, r2, r4, r6, dvlj;
	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
			);
			r2 = r * r;
			r4 = r2 * r2;
			r6 = r4 * r2;
			dvlj = 4.0e0 * (-1 / (r4*r4) + 1 / (r6*r6*r2));
			gradient[i] += dvlj * (x[i] - x[j]);
			gradient[j] -= dvlj * (x[i] - x[j]);
			gradient[i + natm] += dvlj * (x[i + natm] - x[j + natm]);
			gradient[j + natm] -= dvlj *(x[i + natm] - x[j + natm]);
			gradient[i + 2 * natm] += dvlj * (x[i + 2 * natm] - x[j + 2 * natm]);
			gradient[j + 2 * natm] -= dvlj * (x[i + 2 * natm] - x[j + 2 * natm]);
			
		}
	}
	return true;
}

double Fitness::CoulombEnergy(vector<double> &x, vector<double> & atomsCharge)
{
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	double r;
	double vCoulomb = 0.0e0;
	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
			);
			vCoulomb += (atomsCharge[i] * atomsCharge[j]) / r;
		}
	}
	if (isnan(vCoulomb))
		return 1.0e99;
	return vCoulomb;
}


void Fitness::CoulombGradient(
	vector<double> &x, 
	vector<double> &gradient,
	vector<double> &atomsCharge)
{
	// x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3
	int natm = x.size() / 3;
	for (size_t i = 0; i < gradient.size(); i++)
		gradient[i] = 0.0e0;

	double r, dCoulomb;
	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
			);
			dCoulomb = (atomsCharge[i] * atomsCharge[j]) / (r * r * r);
			gradient[i] += dCoulomb * (x[i] - x[j]);
			gradient[j] -= dCoulomb * (x[i] - x[j]);
			gradient[i + natm] += dCoulomb * (x[i + natm] - x[j + natm]);
			gradient[j + natm] -= dCoulomb *(x[i + natm] - x[j + natm]);
			gradient[i + 2 * natm] += dCoulomb * (x[i + 2 * natm] - x[j + 2 * natm]);
			gradient[j + 2 * natm] -= dCoulomb * (x[i + 2 * natm] - x[j + 2 * natm]);
		}
	}
}

void Fitness::CoulombGradient(
	const state_type &x,
	state_type &dxdt)
{
	/* state_type:
	xe1 - 0, xp1 - 1, xe2 - 2, xp2 - 3,
	ye1 - 4, yp1 - 5, ye2 - 6, yp2 - 7,
	ze1 - 8, zp1 - 9, ze2 - 10, zp2 - 11
	vxe1 - 12, vxp1 - 13, vxe2 - 14, vxp2 - 15,
	vye1 - 16, vyp1 - 17, vye2 - 18, vyp2 - 19,
	vze1 - 20, vzp1 - 21, vze2 - 22, vzp2 - 23 */

	int natm = x.size() / 6;
	for (size_t i = 0; i < dxdt.size(); i++)
		dxdt[i] = 0.0e0;

	for (size_t i = 0; i < 12; i++)
		dxdt[i] = x[i + 12]; // velocities dxdt-x1,x2 stays at x12,x13...

	//dxdt de v ==> F/m
	double r, dCoulomb;
	for (int i = 0; i < (natm - 1); i++)
	{
		for (int j = (i + 1); j < natm; j++)
		{
			r = sqrt(
				(x[i] - x[j])*(x[i] - x[j]) +
				(x[i + natm] - x[j + natm])*(x[i + natm] - x[j + natm]) +
				(x[i + 2 * natm] - x[j + 2 * natm])*(x[i + 2 * natm] - x[j + 2 * natm])
			);
			dCoulomb = (atomsChargeFit[i] * atomsChargeFit[j]) / (r * r * r);
			dxdt[12 + i] += dCoulomb * (x[i] - x[j]);
			dxdt[12 + j] -= dCoulomb * (x[i] - x[j]);
			dxdt[12 + i + natm] += dCoulomb * (x[i + natm] - x[j + natm]);
			dxdt[12 + j + natm] -= dCoulomb * (x[i + natm] - x[j + natm]);
			dxdt[12 + i + 2 * natm] += dCoulomb * (x[i + 2 * natm] - x[j + 2 * natm]);
			dxdt[12 + j + 2 * natm] -= dCoulomb * (x[i + 2 * natm] - x[j + 2 * natm]);
		}
	}

	dxdt[13] /= mProton;
	dxdt[15] /= mProton;
	dxdt[17] /= mProton;
	dxdt[19] /= mProton;
	dxdt[21] /= mProton;
	dxdt[23] /= mProton;

}

void Fitness::calculateTotalCoulombSystemEnergy(
	std::vector<double> &x,
	std::vector<double> &v,
	std::vector<double> &atomsCharge,
	std::vector<double> &atomsMass,
	std::ofstream &printFile_)
{
	double ePotential = CoulombEnergy(x, atomsCharge);
	int natm = atomsCharge.size();
	double eKinetic = 0.0e0;
	for (int i = 0; i < natm; i++)
	{
		eKinetic += 0.5e0 * atomsMass[i] * (
			(v[i] * v[i]) +
			(v[i + natm] * v[i + natm]) +
			(v[i + 2 * natm] * v[i + 2 * natm]));
	}
	double eTot = ePotential + eKinetic;
	printFile_ << fixed << setprecision(12) << setw(16) << ePotential << ";"
		<< fixed << setprecision(12) << setw(16) << eKinetic << ";"
		<< fixed << setprecision(12) << setw(16) << eTot << ";";
}


double Fitness::calculateTotalEnergyCoulomb(
	std::vector<double> &x,
	std::vector<double> &v,
	std::vector<double> &atomsCharge,
	std::vector<double> &atomsMass)
{
	double ePotential = CoulombEnergy(x, atomsCharge);
	int natm = atomsCharge.size();
	double eKinetic = 0.0e0;
	for (int i = 0; i < natm; i++)
	{
		eKinetic += 0.5e0 * atomsMass[i] * (
			(v[i] * v[i]) +
			(v[i + natm] * v[i + natm]) +
			(v[i + 2 * natm] * v[i + 2 * natm]));
	}
	return ePotential + eKinetic;
}

double Fitness::calculateTotalAngularMomentum(
	std::vector<double> &x,
	std::vector<double> &v,
	std::vector<double> &atomsMass)
{
	/* x and v
	xe1 - 0, xp1 - 1, xe2 - 2, xp2 - 3,
		ye1 - 4, yp1 - 5, ye2 - 6, yp2 - 7,
		ze1 - 8, zp1 - 9, ze2 - 10, zp2 - 11
	*/

	AuxMath auxMath_;
	vector<double> Le1, Le2, Lp1, Lp2, Ltotal;

	Le1 = auxMath_.vectorProduct(
		x[0], x[4], x[8],
		v[0], v[4], v[8]);
	Lp1 = auxMath_.vectorProduct(
		x[1], x[5], x[9],
		v[1], v[5], v[9]);
	Le2 = auxMath_.vectorProduct(
		x[2], x[6], x[10],
		v[2], v[6], v[10]);
	Lp2 = auxMath_.vectorProduct(
		x[3], x[7], x[11],
		v[3], v[7], v[11]);

	Ltotal.resize(3);
	Ltotal[0] = Le1[0] + mProton * Lp1[0] + Le2[0] + mProton * Lp2[0];
	Ltotal[1] = Le1[1] + mProton * Lp1[1] + Le2[1] + mProton * Lp2[1];
	Ltotal[2] = Le1[2] + mProton * Lp1[2] + Le2[2] + mProton * Lp2[2];

	return auxMath_.norm(Ltotal[0], Ltotal[1], Ltotal[2]);
}


void Fitness::printCenterOfMass(
	std::vector<double> &x,
	std::vector<double> atomsMass,
	std::ofstream &printFile_)
{
	double cmx = 0.0e0;
	double cmy = 0.0e0;
	double cmz = 0.0e0;

	int natm = x.size() / 3;
	for (int i = 0; i < natm; i++)
	{
		cmx += x[i] * atomsMass[i];
		cmy += x[i + natm] * atomsMass[i + natm];
		cmz += x[i + 2 * natm] * atomsMass[i + 2 * natm];
	}
	cmx /= natm;
	cmy /= natm;
	cmz /= natm;

	printFile_ << fixed << setprecision(12) << setw(16) << cmx << ";"
		<< fixed << setprecision(12) << setw(16) << cmy << ";"
		<< fixed << setprecision(12) << setw(16) << cmz << ";";
}


double Fitness::optimizeLennardJones(std::vector<double> &x, int fitType)
{
#ifdef useDlib
	using namespace dlib;

	int size = x.size();
	column_vector starting_point(size);
	for (int i = 0; i < size; i++)
		starting_point(i) = x[i];

	double fMin = find_min(bfgs_search_strategy(),
		objective_delta_stop_strategy(1e-6),
		FunctionDlib(size, fitType),
		DerivativeDlib(size, fitType),
		starting_point,
		-1.0e99);

	for (int i = 0; i < size; i++)
		x[i] = starting_point(i);

	return fMin;
#else
	return lennardJones(x);
#endif	
}

/* EXMPLO DE OPTIMIZE
main:
InitializeAtoms init_;
vector<double> x = init_.generateCluster(20, 0.2, 2.5);
Fitness fit_;
printAtomsVectorDouble(x, "teste1.xyz");
fit_.optimizeLennardJones(x, 0);
printAtomsVectorDouble(x, "teste2.xyz");
*/
