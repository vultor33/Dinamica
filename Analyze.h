#ifndef ANALYZE_H
#define ANALYZE_H

#include <vector>
#include <string>

#include "Coordstructs.h"

class Analyze
{
public:
	Analyze();
	~Analyze();

	void takeChargeDistribution();

	void takeReactedResults();

private:
	std::vector< std::vector<CoordXYZ> > readSimulationInput(std::string fileName);

	void addToCount(std::vector<int> & countHisto, double proj, double rpp);

	std::string takeStringUntilCharacter(std::string &entryString, std::string referenceCharacter);

	void printCoulombAtoms(std::vector<double> & atoms, std::string testName, std::vector<double> &atomsCharge);

	void chargeDistribution(std::string fileName);

	//data
	double bondDistanceCut;




};


#endif

