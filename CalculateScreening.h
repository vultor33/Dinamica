#ifndef CALCULATESCREENING_H
#define CALCULATESCREENING_H

class CalculateScreening
{
public:
	CalculateScreening();

	~CalculateScreening();

	/*
	void calcBohr();
	void calcBohrEllipse();
	void calcBohrEllipseAngle();
	void calcEllipseLh2b();
	*/

	void screenDynamicCenter(int symmetricType, int initialPositions, int kOnly = -1);

	void screenDynamicRear(int symmetricType, int initialPositions, int kOnly = -1);

	void calcOne();



};


#endif