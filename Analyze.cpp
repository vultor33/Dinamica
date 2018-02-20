#include "Analyze.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Coordstructs.h"
#include "AuxMath.h"

using namespace std;

Analyze::Analyze()
{
	bondDistanceCut = 2.5e0;
}

Analyze::~Analyze(){}


void Analyze::takeChargeDistribution()
{
	system("ls -1 > allNames.txt");
	ifstream names_("allNames.txt");
	string name;

	while (getline(names_, name))
	{
		string simulFlag = takeStringUntilCharacter(name, ".");
		if (name == "xyz")
		{
			string fullName = simulFlag + "." + name;
			chargeDistribution(fullName);
		}
	}


}


void Analyze::chargeDistribution(string fileName)
{
	vector< vector<CoordXYZ> > allMol = readSimulationInput(fileName);
	ofstream results_;
	//results_.open("electronPlot.csv", std::ofstream::out | std::ofstream::app);
	results_.open("electronPlot.csv");
	ofstream rppGraph_("rppGraph.csv");
	ofstream yTraject_("yTrajectory.csv");
	vector<int> countHisto(12);
	for (size_t i = 0; i < countHisto.size(); i++)
		countHisto[i] = 0;
	double rppMean = 0.0e0;
	double rpp0 = 0.0e0;
	vector<double> initialEP1;
	int kRpp = 1;
	double oldAmpl = -1.0e0;
	int iOldAmpl = -1;
	double midAmpl = -1.0e0;
	int iMidAmpl = -1;
	vector<int> amplMax;
	for (size_t i = 0; i < allMol.size(); i++)
	{
		double rpp = sqrt(
			(allMol[i][1].x - allMol[i][3].x)*(allMol[i][1].x - allMol[i][3].x)
			+ (allMol[i][1].y - allMol[i][3].y)*(allMol[i][1].y - allMol[i][3].y)
			+ (allMol[i][1].z - allMol[i][3].z)*(allMol[i][1].z - allMol[i][3].z));
		rppMean += rpp;

		/* plot rpp all
		if (i % 100 == 0)
		{
			rppGraph_ << kRpp << " ; " << rpp << endl;
			kRpp++;
		}
		*/

//		if (rpp < bondDistanceCut) //apagar

		vector<double> vecPP(3);
		vecPP[0] = -allMol[i][3].x + allMol[i][1].x;
		vecPP[1] = -allMol[i][3].y + allMol[i][1].y;
		vecPP[2] = -allMol[i][3].z + allMol[i][1].z;
		vector<double> vecEP1(3);
		vecEP1[0] = -allMol[i][0].x + allMol[i][1].x;
		vecEP1[1] = -allMol[i][0].y + allMol[i][1].y;
		vecEP1[2] = -allMol[i][0].z + allMol[i][1].z;
		vector<double> vecEP2(3);
		vecEP2[0] = -allMol[i][2].x + allMol[i][1].x;
		vecEP2[1] = -allMol[i][2].y + allMol[i][1].y;
		vecEP2[2] = -allMol[i][2].z + allMol[i][1].z;


		if (i == 0)
		{
			rpp0 = rpp;
		}


		AuxMath auxMath_;
		double prodEsc = auxMath_.escalarProduct(
			vecPP[0], vecPP[1], vecPP[2],
			vecEP1[0], vecEP1[1], vecEP1[2]);
		double prodEsc2 = auxMath_.escalarProduct(
			vecPP[0], vecPP[1], vecPP[2],
			vecEP2[0], vecEP2[1], vecEP2[2]);

		addToCount(countHisto, prodEsc / rpp, rpp);
		addToCount(countHisto, prodEsc2 / rpp, rpp);

		// y projection
		// proj vector


		double projJOnEpp = vecPP[1] / rpp;
		vector<double> vProjJOnEpp(3);
		vProjJOnEpp[0] = projJOnEpp * vecPP[0] / rpp;
		vProjJOnEpp[1] = projJOnEpp * vecPP[1] / rpp;
		vProjJOnEpp[2] = projJOnEpp * vecPP[2] / rpp;
		vector<double> ppPerpendicular(3);
		ppPerpendicular[0] = -vProjJOnEpp[0];
		ppPerpendicular[1] = 1.0e0 - vProjJOnEpp[1];
		ppPerpendicular[2] = -vProjJOnEpp[2];
		double normPerp = auxMath_.norm(ppPerpendicular[0], ppPerpendicular[1], ppPerpendicular[2]);

		double yProj = auxMath_.escalarProduct(
			ppPerpendicular[0], ppPerpendicular[1], ppPerpendicular[2],
			vecEP1[0], vecEP1[1], vecEP1[2]) / normPerp;

		yTraject_ << ((prodEsc / rpp) - rpp0 / 2.0e0) << "  " << yProj << endl;
		rppGraph_ << kRpp << "  " << rpp << endl;
		kRpp++;

		if ((midAmpl > oldAmpl) && (midAmpl > rpp))
		{
			//maximo
			amplMax.push_back(iMidAmpl);
		}
		oldAmpl = midAmpl;
		iOldAmpl = iMidAmpl;
		midAmpl = rpp;
		iMidAmpl = kRpp;
	}

	results_ << "-0.5-  " << countHisto[0] << endl;
	results_ << "-0.5a-0.3  " << countHisto[1] << endl;
	results_ << "-0.3a-0.1  " << countHisto[2] << endl;
	results_ << "-0.1a0.1  " << countHisto[3] << endl;
	results_ << "0.1a0.3  " << countHisto[4] << endl;
	results_ << "0.3a0.5  " << countHisto[5] << endl;
	results_ << "0.5a0.7  " << countHisto[6] << endl;
	results_ << "0.7a0.9  " << countHisto[7] << endl;
	results_ << "0.9a1.1  " << countHisto[8] << endl;
	results_ << "1.1a1.3  " << countHisto[9] << endl;
	results_ << "1.3a1.5  " << countHisto[10] << endl;
	results_ << "1.5+  " << countHisto[11] << endl;
	results_ << endl;
	rppGraph_.close();

	cout << endl << endl;
	double amplMeanValue = 0.0e0;
	for (size_t i = 1; i < amplMax.size(); i++)
	{
		amplMeanValue += (amplMax[i] - amplMax[i - 1]);
	}
	amplMeanValue /= (double)(amplMax.size() - 1);

	results_ << "r mean: " << rppMean / allMol.size() << endl;
	results_ << "periodo oscilacao:  " << amplMeanValue << endl;




	
	
	
	
	results_.close();
	





}

vector< vector<CoordXYZ> > Analyze::readSimulationInput(string fileName)
{
	ifstream file_(fileName.c_str());
	string line;
	vector< vector<CoordXYZ> > allMol;
	while (getline(file_, line))
	{
		stringstream auxConv;
		auxConv << line;
		int nAtoms;
		auxConv >> nAtoms;
		getline(file_, line);
		vector<CoordXYZ> mol;
		for (int i = 0; i < nAtoms; i++)
		{
			getline(file_, line);
			stringstream auxConv2;
			auxConv2 << line;
			CoordXYZ auxCoord;
			auxConv2 >> auxCoord.atomlabel
				>> auxCoord.x
				>> auxCoord.y
				>> auxCoord.z;
			mol.push_back(auxCoord);
		}
		allMol.push_back(mol);
	}
	return allMol;
}


void Analyze::addToCount(vector<int> &countHisto, double proj, double rpp)
{
	if (proj < -0.5e0 * rpp)
		countHisto[0]++;
	else if ((proj > -0.5e0 * rpp) && (proj < -0.3e0 * rpp))
		countHisto[1]++;
	else if ((proj > -0.3e0 * rpp) && (proj < -0.1e0 * rpp))
		countHisto[2]++;
	else if ((proj > -0.1e0 * rpp) && (proj < 0.1e0 * rpp))
		countHisto[3]++;
	else if ((proj > 0.1e0 * rpp) && (proj < 0.3e0 * rpp))
		countHisto[4]++;
	else if ((proj > 0.3e0 * rpp) && (proj < 0.5e0 * rpp))
		countHisto[5]++;
	else if ((proj > 0.5e0 * rpp) && (proj < 0.7e0 * rpp))
		countHisto[6]++;
	else if ((proj > 0.7e0 * rpp) && (proj < 0.9e0 * rpp))
		countHisto[7]++;
	else if ((proj > 0.9e0 * rpp) && (proj < 1.1e0 * rpp))
		countHisto[8]++;
	else if ((proj > 1.1e0 * rpp) && (proj < 1.3e0 * rpp))
		countHisto[9]++;
	else if ((proj > 1.3e0 * rpp) && (proj < 1.5e0 * rpp))
		countHisto[10]++;
	else if (proj > 1.5e0 * rpp)
		countHisto[11]++;

}






void Analyze::takeReactedResults()
{
	system("ls -1 > allNames.txt");
	ifstream names_("allNames.txt");
	ofstream movieSimul_("movie-simul.csv");
	movieSimul_ << "seed ; Temperature ; Impact Parameter ; Method ; Proton Distance" << endl;
	string name;
	double rProtonsCut = 2.5e0;
	double eCut = 50.0e0;
	vector<double> atomsCharge(4);
	atomsCharge[0] = -1.0e0;
	atomsCharge[1] = 1.0e0;
	atomsCharge[2] = -1.0e0;
	atomsCharge[3] = 1.0e0;

	while (getline(names_, name))
	{
		string simulFlag = takeStringUntilCharacter(name, "-");
		if (simulFlag == "simulation")
		{
			string fullName = "simulation-" + name;

			// seed, T, b, metodo
			string seedString = takeStringUntilCharacter(name, "-");
			string tempString = takeStringUntilCharacter(name, "-");
			string impactString = takeStringUntilCharacter(name, "-");
			string methodString = takeStringUntilCharacter(name, ".");

			ifstream enterSimul_(fullName.c_str());
			string simulLine;
			getline(enterSimul_, simulLine);
			vector<double> x(12);
			vector< vector<double> > allX;
			while (getline(enterSimul_, simulLine))
			{
				stringstream convert;
				convert << simulLine;
				string dummy;
				for (int i = 0; i < 12; i++)
				{
					convert >> x[i];
					convert >> dummy;
					convert >> dummy;
					convert >> dummy;
				}
				allX.push_back(x);
			}
			double rProtons = sqrt(
				(x[1] - x[3])*(x[1] - x[3])
				+ (x[5] - x[7])*(x[5] - x[7])
				+ (x[9] - x[11])*(x[9] - x[11]));

			double rE1 = sqrt(
				(x[0] * x[0])
				+ (x[4] * x[4])
				+ (x[8] * x[8]));

			double rE2 = sqrt(
				(x[2] * x[2])
				+ (x[6] * x[6])
				+ (x[10] * x[10]));

			if ((rProtons < rProtonsCut) && (rE1 < eCut) && (rE2 < eCut))
			{
				for (size_t i = 0; i < allX.size(); i++)
				{
					printCoulombAtoms(allX[i], methodString + "-" + seedString + "-" + tempString + "-" + impactString + ".xyz", atomsCharge);
				}

				movieSimul_ << seedString << " ; "
					<< tempString << " ; "
					<< impactString << " ; "
					<< methodString << " ; "
					<< rProtons << endl;
			}
		}
	}

}




string Analyze::takeStringUntilCharacter(string &entryString, string referenceCharacter)
{
	size_t hyfPos = entryString.find(referenceCharacter);
	string before = entryString.substr(0, hyfPos);
	entryString = entryString.substr(hyfPos + 1, entryString.size());
	return before;
}


void Analyze::printCoulombAtoms(vector<double> & atoms, string testName, vector<double> &atomsCharge)
{
	int natm = atoms.size() / 3;
	ofstream teste_;
	teste_.open(testName.c_str(), std::ofstream::out | std::ofstream::app);
	teste_ << natm << endl << "t" << endl;
	for (int i = 0; i < natm; i++)
	{
		if (atomsCharge[i] == 1.0e0)
			teste_ << "Au ";
		else
			teste_ << "H ";

		teste_ << atoms[i] << "  "
			<< atoms[i + natm] << "  "
			<< atoms[i + 2 * natm] << endl;
	}
	teste_.close();
}





