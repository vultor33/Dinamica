#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>
#include <fstream>

class Fitness
{
public:
	Fitness();
	~Fitness();

	double fit(std::vector<double> &point, int type);

	double optimizeLennardJones(std::vector<double> &x, int fitType);

	double runGamess(
		std::vector<double> &x,
		std::vector< std::string > &options,
		std::string gamessPath,
		std::string gamessScr,
		std::string nProc);

	bool lennardJonesGradient(std::vector<double> &x, std::vector<double> &gradient);

	double CoulombEnergy(std::vector<double> &x, std::vector<double> & atomsCharge);

	bool CoulombGradient(
		std::vector<double> & x,
		std::vector<double> & gradient,
		std::vector<double> &atomsCharge);

	void calculateTotalCoulombSystemEnergy(
		std::vector<double> &x,
		std::vector<double> &v,
		std::vector<double> &atomsCharge,
		std::vector<double> &atomsMass,
		std::ofstream &printFile_);

	void printCenterOfMass(
		std::vector<double> &x, 
		std::vector<double> atomsMass,
		std::ofstream &printFile_);

private:
	double lennardJones(std::vector<double> &x);

};

#endif
