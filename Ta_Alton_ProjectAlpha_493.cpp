// ConsoleApplication14.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <assert.h>
#include <string.h>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>

using namespace std;

#define ATRAND (double)rand()/RAND_MAX

double generateGaussianNoise(double mu, double sigma)																					// taken from wikipedia box muller
{
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

double average(vector<double>* pv) {
	double sum = 0;
	for (int i = 0; i < pv->size(); i++) {
		sum = sum + pv->at(i);
	}
	return sum / pv->size();
};

double stdev() {

};

class arm {
public:
	double reward = 0;
	double mu = 0;
	double sigma = 0;
	double V_t = 0;
	double alpha = 0;
	void init();
	void pull();
	void update();
};

void arm::init() {
	double mu = rand() % 100;
	double sigma = rand() % 5;
};

void arm::pull() {
	reward = generateGaussianNoise(mu, sigma)
}

void arm::update() {
	V_t = reward*alpha + V_t*(1 - alpha);
}

class agent {
public:																																	// learner
	double alpha = 0;
	double epsilon = 0;
	int action = 0;
	int V_t_1 = 0;
	int V_t_2 = 0;
	int V_t_3 = 0;
	int arm_pulled = 0;
	void init();
	void decide();
	void act();
	void react();
};

void agent::init() {
	double alpha = 0.1;																													// set as such because it's deemed as a good value for alpha
	double epsilon = 0.7;																												// set as such to encourage exploration at the start
};

void agent::decide() {
	if (epsilon >= ATRAND)																												// compares epsilon to a random number between [0,1]
		action = 1;																											// if greater or equal, random action taken
	else
		action = 2;																											// if lesser, greedy option taken
}

void agent::act() {
	if (action == 1) {																											// pull a random arm, exploration
		int decide = rand() % 3 + 1;
		if (decide == 1) {																												// if random number = 1, first arm pulled
			pull.arm_a;
			arm_pulled = 1;
		}
		if (decide == 2) {																												// if random number = 2, second arm pulled
			pull.arm_b;
			arm_pulled = 2;
		}
		if (decide == 3) {																												// if random number = 3, third arm pulled
			pull.arm_c;
			arm_pulled = 3;
		}
		else
			assert(1 == 0);																												// something went wrong!
	}

	if (action == 2) {																																	// pick greatest value
		if ((V_t_1 > V_t_2) && (V_t_1 > V_t_3)) {																		// if arm 1 > 2 and 3
			pull.arm_a;
			arm_pulled = 1;
		}
		if ((V_t_2 > V_t_1) && (V_t_2 > V_t_3)) {																		// if arm 2 > 1 and 3
			pull.arm_b;
			arm_pulled = 2;
		}
		if ((V_t_3 > V_t_1) && (V_t_3 > V_t_1)) {																		// if arm 3 > 1 and 2
			pull.arm_c;
			arm_pulled = 3;

		}
		else {																															// for first pull, if epsilon is not random (this is a hard coded action)
			if (decide == 1) {
				pull.arm_a;
				arm_pulled = 1;
			}
			if (decide == 2) {
				pull.arm_b;
				arm_pulled = 2;
			}
			if (decide == 3) {
				pull.arm_c;
				arm_pulled = 3;
			}
			else
				assert(1 == 0);
		}


	}
}


void agent::react() {
	if (arm_pulled = 1) {
		V_t_1 = reward * alpha + V_t_1 * (1 - alpha);
	}
	if (arm_pulled = 2) {
		V_t_2 = reward * alpha + V_t_2 * (1 - alpha);
	}
	if (arm_pulled = 3) {
		V_t_2 = reward * alpha + V_t_2 * (1 - alpha);
	}
}

class MAB {
public:
	void init();
};


void MAB::init() {
	arm arm_a;
	arm arm_b;
	arm arm_c;
	arm_a.init;
	arm_b.init;
	arm_c.init;
};


int main()
{
	srand(time(NULL));
	int test = 0;
	cout << test << endl;
	vector<arm> arms;
	int n = 3;
	for (int i = 0; i < n; i++) {
		arm A;
		A.init;
		arms.push_back(A);
	}																																			// make n arms
	return 0;
}
