#pragma once
#ifndef ODEINTSAFEOBSERVER_H
#define ODEINTSAFEOBSERVER_H

#include <vector>

typedef std::vector< double > state_type;

class CounterWrite
{
public:
	CounterWrite(int maxCounter_in)
		: maxCounter(maxCounter_in)
	{
		counter = 0;
	}

	~CounterWrite() {}

	void addCounter()
	{
		counter++;
		if (counter > maxCounter)
		{
			throw 1;
		}
	}

private:
	int maxCounter;
	int counter;
};

struct write_state
{
	CounterWrite *counter;

	void operator()(const state_type &, const double) const
	{
		counter->addCounter();
	}
};



#endif

