#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

string takeStringUntilCharacter(string &entryString, string referenceCharacter);

void printCoulombAtoms(vector<double> & atoms, string testName, vector<double> &atomsCharge);

int main()
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

	while(getline(names_,name))
	{
		string simulFlag = takeStringUntilCharacter(name,"-");
		if(simulFlag == "simulation")
		{
			string fullName = "simulation-" + name;

			// seed, T, b, metodo
			string seedString = takeStringUntilCharacter(name,"-");
			string tempString = takeStringUntilCharacter(name,"-");
			string impactString = takeStringUntilCharacter(name,"-");
			string methodString = takeStringUntilCharacter(name,".");
			
			ifstream enterSimul_(fullName.c_str());
			string simulLine;
			getline(enterSimul_,simulLine);
			vector<double> x(12);
			vector< vector<double> > allX;
			while(getline(enterSimul_,simulLine))
			{
				stringstream convert;
				convert << simulLine;
				string dummy;
				for(int i = 0; i < 12; i++)
				{
					convert >> x[i];
					convert >> dummy;
					convert >> dummy;
					convert >> dummy;
				}
				allX.push_back(x);
			}
			double rProtons = sqrt(
				(x[1]-x[3])*(x[1]-x[3])
				+(x[5]-x[7])*(x[5]-x[7])
				+(x[9]-x[11])*(x[9]-x[11]));

			double rE1 = sqrt(
				(x[0]*x[0])
				+(x[4]*x[4])
				+(x[8]*x[8]));

			double rE2 = sqrt(
				(x[2]*x[2])
				+(x[6]*x[6])
				+(x[10]*x[10]));

			if((rProtons < rProtonsCut)&&(rE1 < eCut)&&(rE2 < eCut))
			{
				for(size_t i = 0; i < allX.size(); i++)
				{
					printCoulombAtoms(allX[i],methodString + "-" + seedString + "-" + tempString + "-" + impactString + ".xyz",atomsCharge);
				}	

				movieSimul_ << seedString << " ; "
					<< tempString << " ; "
					<< impactString << " ; "
					<< methodString << " ; "
					<< rProtons << endl;
			}				
		}
	}

	return 0;

}



string takeStringUntilCharacter(string &entryString, string referenceCharacter)
{
	size_t hyfPos = entryString.find(referenceCharacter);
	string before = entryString.substr(0,hyfPos);
	entryString = entryString.substr(hyfPos + 1, entryString.size());
	return before;
}


void printCoulombAtoms(vector<double> & atoms, string testName, vector<double> &atomsCharge)
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


