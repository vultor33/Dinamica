#ifndef READDYMINPUT_H
#define READDYMINPUT_H

#include "DynamicsStructs.h"

class ReadDymInput
{
public:
	ReadDymInput();

	~ReadDymInput();

	void defineMethodBohr(double angle);

	void electronPlot();

	void analyzePlot();

	void addIToName(int i);

	DymOptions getDymOptions();

private:
	DymOptions dymOptions_;

	void generateDefaultOptions();

	double initialVelocityKinecticTheory(double tempKelvin);

};


#endif