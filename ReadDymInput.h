#ifndef READDYMINPUT_H
#define READDYMINPUT_H

#include "DynamicsStructs.h"

class ReadDymInput
{
public:
	ReadDymInput();

	~ReadDymInput();

	void defineMethodSymmetries(
		int symm,
		int initialPositionType,
		double rElec,
		double rProton,
		double angle);

	/*
	void defineMethodBohr(double angle);
	void defineMethodEllipseBohr(double rElec, double dProton);
	void defineMethodEllipseBohrAngle(double rElec, double dProton, double angle);
	void defineMethodLh2b(double rElec, double rProton);
	*/

	void electronPlot();

	void analyzePlot();

	void setSeed(int seed);

	void addIToName(int i);

	DymOptions getDymOptions();

private:
	DymOptions dymOptions_;

	void generateDefaultOptions();

	double initialVelocityKinecticTheory(double tempKelvin);

};


#endif