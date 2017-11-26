#ifndef GENERATEATOM_H
#define GENERATEATOM_H

#include <vector>
#include <string>

class GenerateAtom
{
public:
	GenerateAtom(int seed);
	~GenerateAtom();

	void generateTwoRandomAtoms(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge);

	void generateTwoSymmetricAtoms(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge);

	void generateTwoIdenticalSymmetricAtoms(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge);

	void generateTwoAntiSymmetricAtoms(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge);

	void generateTwoIdenticalAntiSymmetricAtoms(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge);

	void translateToCenterOfMass(std::vector<double> &x, std::vector<double> &atomsMass);

	// NOT WORKING
	void velocityCmCorrection(std::vector<double> &v, std::vector<double> &atomsMass);//4 particles SPECIFIC

private:
	double mProton;
	double pi_;

	// Le = 1 ; P = 0
	void generateInitialPositionAndVelocity(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge);

	std::vector<double> addAtom(
		std::vector<double> & property1,
		std::vector<double> & property2);

	std::vector<double> insertVector(
		std::vector<double> & vec1,
		std::vector<double> & vec2);

	std::vector<double> unitarySphericalVector();

	double randomNumber(double fMin, double fMax);

};


#endif