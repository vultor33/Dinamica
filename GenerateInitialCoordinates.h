#ifndef GENERATEINITIALCOORDINATES_H
#define GENERATEINITIALCOORDINATES_H

#include <vector>
#include <string>

#include "DynamicsStructs.h"

class GenerateInitialCoordinates
{
public:
	GenerateInitialCoordinates();

	~GenerateInitialCoordinates();

	bool generateInitial(
		DymOptions &dymOptions_,
		std::vector<double> &x, 
		std::vector<double> &v, 
		std::vector<double> &atomsMass, 
		std::vector<double> &atomsCharge);

private:
	double mProton;
	double pi_;



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

	void generateBohrMolecule(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge,
		double energy,
		double angle
		);

	void generateBohrEllipseMolecule(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge,
		double energy,
		double rElec,
		double rProton);

	void generateBohrEllipseMoleculeAngle(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge,
		double energy,
		double rElec,
		double rProton,
		double angle);

	void calcBohrParams(
		double energy,
		double &rProtonInit,
		double &rElecInit,
		double &vInit);

	void translateToCenterOfMass(std::vector<double> &x, std::vector<double> &atomsMass);



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


	// NOT WORKING
	void velocityCmCorrection(std::vector<double> &v, std::vector<double> &atomsMass);//4 particles SPECIFIC









	bool generateElectronsAtCenter(
		std::vector<double> &xPositions,
		std::vector<double> &vVelocities,
		std::vector<double> &atomsMass,
		std::vector<double> &atomCharge,
		DymOptions &dymOptions);


	void initialVelocities(
		double deltaEnergy,
		std::vector<double> &v,
		DymOptions &dymOption_);





};


#endif