#pragma once
#ifndef FITNESSCOULOMBODEINT_H
#define FITNESSCOULOMBODEINT_H

#include <iostream>

#include "Fitness.h"
#include "TrajectorySymmetrizer.h"

typedef std::vector< double > state_type;
/* state_type:
xe1 - 0, xp1 - 1, xe2 - 2, xp2 - 3,
ye1 - 4, yp1 - 5, ye2 - 6, yp2 - 7,
ze1 - 8, zp1 - 9, ze2 - 10, zp2 - 11
vxe1 - 12, vxp1 - 13, vxe2 - 14, vxp2 - 15,
vye1 - 16, vyp1 - 17, vye2 - 18, vyp2 - 19,
vze1 - 20, vzp1 - 21, vze2 - 22, vzp2 - 23 */

//typedef boost::numeric::odeint::runge_kutta_dopri5< double > stepper_type;
//typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;


class FitnessCoulombOdeint
{
public:

	FitnessCoulombOdeint(int symmetrizeOption_in)
		: symmetrizeOption(symmetrizeOption_in) {}

	void operator() (const state_type &x, state_type &dxdt, const double t)
	{
		// x[0] = x --- dxdt[0] = x[1] velocity x
		// x[1] = p --- dxdt[1] = a = F/m 
		fit_.CoulombGradient(x, dxdt);
		trSymm_.symmetrize(symmetrizeOption, dxdt);
	}

private:
	Fitness fit_;
	TrajectorySymmetrizer trSymm_;

	int symmetrizeOption;

};

#endif
