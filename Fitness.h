#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>
#include <fstream>

typedef std::vector< double > state_type;

class Fitness
{
public:
	Fitness();
	~Fitness();

	double fit(std::vector<double> &point, int type);

	double optimizeLennardJones(std::vector<double> &x, int fitType);

	bool lennardJonesGradient(std::vector<double> &x, std::vector<double> &gradient);

	double CoulombEnergy(std::vector<double> &x, std::vector<double> & atomsCharge);

	void CoulombGradient(
		std::vector<double> & x,
		std::vector<double> & gradient,
		std::vector<double> &atomsCharge);

	void CoulombGradient(
		const state_type &x,
		state_type &dxdt);

	void calculateTotalCoulombSystemEnergy(
		std::vector<double> &x,
		std::vector<double> &v,
		std::vector<double> &atomsCharge,
		std::vector<double> &atomsMass,
		std::ofstream &printFile_);

	double calculateTotalEnergyCoulomb(
		std::vector<double> &x,
		std::vector<double> &v,
		std::vector<double> &atomsCharge,
		std::vector<double> &atomsMass);

	double calculateTotalAngularMomentum(
		std::vector<double> &x,
		std::vector<double> &v,
		std::vector<double> &atomsMass);

	void printCenterOfMass(
		std::vector<double> &x, 
		std::vector<double> atomsMass,
		std::ofstream &printFile_);

private:
	double lennardJones(std::vector<double> &x);

	double mProton;
	std::vector<double> atomsChargeFit;
//	std::vector<double> atomsMassFit;

};

#endif
