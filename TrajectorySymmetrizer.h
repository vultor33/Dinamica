#ifndef TRAJECTORYSYMMETRIZER_H
#define TRAJECTORYSYMMETRIZER_H

#include<vector>

typedef std::vector< double > state_type;

class TrajectorySymmetrizer
{
public:
	TrajectorySymmetrizer();

	~TrajectorySymmetrizer();

	void symmetrize(int option,
		std::vector<double> &x,
		std::vector<double> &v);

	void symmetrize(int option,
		state_type &dxdt);

private:

};

#endif

