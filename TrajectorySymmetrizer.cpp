#include "TrajectorySymmetrizer.h"

#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std;

TrajectorySymmetrizer::TrajectorySymmetrizer(){}

TrajectorySymmetrizer::~TrajectorySymmetrizer(){}

void TrajectorySymmetrizer::symmetrize(int option,
	std::vector<double> &x,
	std::vector<double> &v)
{
	// xe1, xp1, xe2, xp2 - 0-3
	// ye1, yp1, ye2, yp2 - 4-7
	// ze1, zp1, ze2, zp2 - 8-11
	switch (option)
	{
	case 0:
		return;

	case 1: // inversion center
		x[2] = -x[0];
		x[3] = -x[1];
		x[6] = -x[4];
		x[7] = -x[5];
		x[10] = -x[8];
		x[11] = -x[9];
		v[2] = -v[0];
		v[3] = -v[1];
		v[6] = -v[4];
		v[7] = -v[5];
		v[10] = -v[8];
		v[11] = -v[9];
		break;

	case 2: // C2x (x->x ; y->-y ; z->-z)
		x[2] = x[0]; //xele2 -> xele1
		x[3] = x[1]; //xpro2 -> xpro1
		x[6] = -x[4]; //yele2 -> -yele1
		x[7] = -x[5]; //ypro2 -> -ypro1
		x[10] = -x[8]; //zele2 -> -zele1
		x[11] = -x[9]; //zpro2 -> -zpro1

		v[2] = v[0]; //vxele2 -> vxele1
		v[3] = v[1]; //vxpro2 -> vxpro1
		v[6] = -v[4]; //vyele2 -> -vyele1
		v[7] = -v[5]; //vypro2 -> -vypro1
		v[10] = -v[8]; //vzele2 -> -vzele1
		v[11] = -v[9]; //vzpro2 -> -vzpro1
		break;

	case 3: // C2z (x-> -x; y->-y  ze2->ze1)
		x[2] = -x[0]; //xele2 -> -xele1
		x[3] = -x[1]; //xpro2 -> -xpro1
		x[6] = -x[4]; //yele2 -> -yele1
		x[7] = -x[5]; //ypro2 -> -ypro1
		x[10] = x[8]; //zele2 -> zele1

		v[2] = -v[0]; //vxele2 -> -vxele1
		v[3] = -v[1]; //vxpro2 -> -vxpro1
		v[6] = -v[4]; //vyele2 -> -vyele1
		v[7] = -v[5]; //vypro2 -> -vypro1
		v[10] = v[8]; //vzele2 -> vzele1
		break;

		
	case 4:  // Cs (x-> x ; y-> -y ; ze2->ze1)
		x[2] = x[0]; //xele2 -> -xele1
		x[3] = x[1]; //xpro2 -> -xpro1
		x[6] = -x[4]; //yele2 -> -yele1
		x[7] = -x[5]; //ypro2 -> -ypro1
		x[10] = x[8]; //zele2 -> zele1

		v[2] = v[0]; //vxele2 -> -vxele1
		v[3] = v[1]; //vxpro2 -> -vxpro1
		v[6] = -v[4]; //vyele2 -> -vyele1
		v[7] = -v[5]; //vypro2 -> -vypro1
		v[10] = v[8]; //vzele2 -> vzele1
		break;


	default:
		cout << "symmetrization not found - exiting" << endl;
		exit(1);
	}

}





