#ifndef TRAJECTORYSYMMETRIZER_H
#define TRAJECTORYSYMMETRIZER_H

#include<vector>

class TrajectorySymmetrizer
{
public:
	TrajectorySymmetrizer();

	~TrajectorySymmetrizer();

	void symmetrize(int option,
		std::vector<double> &x,
		std::vector<double> &v);

private:

};

#endif

