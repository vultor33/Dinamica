#ifndef READDYMINPUT_H
#define READDYMINPUT_H

#include "DynamicsStructs.h"

class ReadDymInput
{
public:
	ReadDymInput();

	~ReadDymInput();

	DymOptions generateDefaultOptions();

private:

	double initialVelocityKinecticTheory(double tempKelvin);

};


#endif