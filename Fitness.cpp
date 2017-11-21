//#define useDlib

#include "Fitness.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include "WriteQuantumInput.h"
#include "ReadQuantumOutput.h"

#ifdef useDlib
#include <dlib/optimization.h>
#include "FunctionDlib.h"
#include "DerivativeDlib.h"
#endif

using namespace std;

Fitness::Fitness(){}

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


bool Fitness::CoulombGradient(
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
	return true;
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
	printFile_ << "pot: " << fixed << setprecision(12) << setw(16) << ePotential
		<< " kin: " << fixed << setprecision(12) << setw(16) << eKinetic
		<< " tot: " << fixed << setprecision(12) << setw(16) << eTot;
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

	printFile_ << " cmx: " << fixed << setprecision(12) << setw(16) << cmx
		<< " cmy: " << fixed << setprecision(12) << setw(16) << cmy
		<< " cmz: " << fixed << setprecision(12) << setw(16) << cmz;
}


double Fitness::runGamess(
	vector<double> &x, 
	vector<string> &options, 
	string gamessPath,
	string gamessScr,
	string nProc)
{
	WriteQuantumInput writeInp_(options);
	int nAtoms = x.size() / 3;
	vector<CoordXYZ> mol(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		mol[i].atomlabel = "N";
		mol[i].x = x[i];
		mol[i].y = x[i + nAtoms];
		mol[i].z = x[i + 2 * nAtoms];
	}
	writeInp_.createInput(mol);

	system(("rm /scr/" + options[1] + "*").c_str());

	system((gamessPath + "  " + options[1] + ".inp  00  " + nProc + " > " + options[1] + ".log").c_str());

	

	ReadQuantumOutput readQ_("gamess");

	readQ_.readOutput((options[1] + ".log").c_str());

	mol = readQ_.getCoordinates();
	for (int i = 0; i < nAtoms; i++)
	{
		x[i] = mol[i].x;
		x[i + nAtoms] = mol[i].y;
		x[i + 2 * nAtoms] = mol[i].z;
	}

	return readQ_.getEnergy();
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
