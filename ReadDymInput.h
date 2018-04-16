#ifndef READDYMINPUT_H
#define READDYMINPUT_H

#include "DynamicsStructs.h"

class ReadDymInput
{
public:
	ReadDymInput();

	~ReadDymInput();

	DymOptions getDymOptions();

	void addIToName(int i);

private:
	DymOptions dymOptions_;

	void generateDefaultOptions();

	double initialVelocityKinecticTheory(double tempKelvin);

};


#endif