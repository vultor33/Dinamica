#ifndef CALCULATESCREENING_H
#define CALCULATESCREENING_H

class CalculateScreening
{
public:
	CalculateScreening();

	~CalculateScreening();

	void screenDynamicCenter(int symmetricType, int initialPositions, int kOnly = -1);

	void screenDynamicCenterPure(int symmetricType, int initialPositions, double iAngle, int kOnly = -1);

	void screenDynamicRear(int symmetricType, int initialPositions, int kOnly = -1);

	void screenDynamicRearPure(int symmetricType, int initialPositions, double iAngle, int kOnly = -1);

	void calcOne();

	void analyzeAll();



};


#endif