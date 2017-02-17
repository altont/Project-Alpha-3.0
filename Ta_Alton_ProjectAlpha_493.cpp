// ConsoleApplication14.cpp : Defines the entry point for the console application.

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
#include <ostream>
#include <iterator>

#define ATRAND (double)rand()/RAND_MAX

using namespace std;

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

double average(vector<double>* pv) {																									// taken from http://www.cplusplus.com/forum/general/42032/
	double sum = 0;
	for (int i = 0; i < pv->size(); i++) {
		sum = sum + pv->at(i);
	}
	return sum / pv->size();
};

double stdev(vector<double>* pv, double average) {
	double E = 0;
	double inv_n = 1.0 / pv->size();
	for (int i = 0; i < pv->size(); i++) {
		E = E + pow(pv->at(i) - average, 2);
	}
	return sqrt(inv_n*E);
};

void testA() {
	int mu = 5;
	int sigma = 1.5;
	double test;
	double avg;
	double variance;

	vector<double> vec;

	for (int n = 0; n < 500; n++) {
		test = generateGaussianNoise(mu, sigma);
		vec.push_back(test);

		avg = average(&vec);

		variance = stdev(&vec, avg);

	}
	cout << "mu is " << mu << endl;
	cout << "avg = " << avg << endl;
	cout << "variance is " << sigma << endl;
	cout << "var = " << variance << endl;
}

class agent {
public:
	void init();
	double epsilon;
	double alpha;
};

void agent::init() {
	alpha = 0.1;
	epsilon = 0.7;
}


class arm {
public:
	void init();
	double mu;
	double sigma;
	double V_t;
};

void arm::init() {
	mu = rand() % 1000 + ATRAND;																											// mean
	sigma = rand() % 50 + ATRAND;																											// variance
	V_t = 0;																															// expected value
}

double pull(arm* parm) {
	double reward;
	reward = generateGaussianNoise(parm->mu, parm->sigma);
	return reward;
}

int decide(agent* plearner) {
	int action;
	double decision = ATRAND;
	if (1 - plearner->epsilon <= decision) {																							// compares if epsilon is less than or equal to a random number
		action = 1;																														// random action
	}
	else {
		action = 2;																														// greedy action
	}
	return action;
}

double act1(int action, double V_t_1, double V_t_2, double V_t_3, arm* parm1, arm* parm2, arm* parm3) {
	int decide;
	int arm_pulled;
	if (action == 1) {
		decide = rand() % 3 + 1;
		if (decide == 1) {
			arm_pulled = 1;
			return arm_pulled;
		}
		if (decide == 2) {
			arm_pulled = 2;
			return arm_pulled;
		}
		if (decide == 3) {
			arm_pulled = 3;
			return arm_pulled;
		}
	}
	if (action == 2) {
		if ((V_t_1 > V_t_2) && (V_t_1 > V_t_3)) {
			arm_pulled = 1;
			return arm_pulled;
		}
		if ((V_t_2 > V_t_1) && (V_t_2 > V_t_3)) {
			arm_pulled = 2;
			return arm_pulled;
		}
		if ((V_t_3 > V_t_1) && (V_t_3 > V_t_1)) {
			arm_pulled = 3;
			return arm_pulled;
		}
		else {
			decide = rand() % 3 + 1;
			if (decide == 1) {
				arm_pulled = 1;
				return arm_pulled;
			}
			if (decide == 2) {
				arm_pulled = 2;
				return arm_pulled;
			}
			if (decide == 3) {
				arm_pulled = 3;
				return arm_pulled;
			}
		}
	}
}

double act2(int action, double V_t_1, double V_t_2, double V_t_3, arm* parm1, arm* parm2, arm* parm3, double arm_pulled) {
	double decide;
	decide = arm_pulled;
	double reward1;
	double reward2;
	double reward3;
	if (action == 1) {
		if (decide == 1) {
			reward1 = pull(parm1);
			return reward1;
		}
		if (decide == 2) {
			reward2 = pull(parm2);
			return reward2;
		}
		if (decide == 3) {
			reward3 = pull(parm3);
			return reward3;
		}
	}
	if (action == 2) {
		if ((V_t_1 > V_t_2) && (V_t_1 > V_t_3)) {
			reward1 = pull(parm1);
			return reward1;
		}
		if ((V_t_2 > V_t_1) && (V_t_2 > V_t_3)) {
			reward2 = pull(parm2);
			return reward2;
		}
		if ((V_t_3 > V_t_1) && (V_t_3 > V_t_2)) {
			reward3 = pull(parm3);
			return reward3;
		}
		else {
			if (decide == 1) {
				reward1 = pull(parm1);
				return reward1;
			}
			if (decide == 2) {
				reward2 = pull(parm2);
				return reward2;
			}
			if (decide == 3) {
				reward3 = pull(parm3);
				return reward3;
			}
		}
	}
}

double react1(double reward1, agent* plearner, arm* parm_1, double alpha, double V_t_1) {
	parm_1->V_t = reward1 * plearner->alpha + V_t_1 * (1 - plearner->alpha);
	return parm_1->V_t;
}

double react2(double reward2, agent* plearner, arm* parm_2, double alpha, double V_t_2) {
	parm_2->V_t = reward2 * plearner->alpha + V_t_2 * (1 - plearner->alpha);
	return parm_2->V_t;
}

double react3(double reward3, agent* plearner, arm* parm_3, double alpha, double V_t_3) {
	parm_3->V_t = reward3 * plearner->alpha + V_t_3 * (1 - plearner->alpha);
	return parm_3->V_t;
}

double update_epsilon(agent* plearner, int i) {
	double e = 2.718181818;
	plearner->epsilon = plearner->epsilon * pow(e, -0.00000099 * i);																// decaying epsilon will allow learner to become more greedy as it learns
	return plearner->epsilon;
}

double update_alpha(agent* plearner, int i) {																					// learning rate increases over pulls
	double e = 2.718181818;
	plearner->alpha = plearner->alpha * pow(e, 0.0000005 * i);																	// should make the learner become greedy faster as
	return plearner->alpha;																										// the growth of alpha > decay of epsilon
}

void testB(arm* parm_ab, arm* parm_bb, arm* parm_cb) {
	parm_ab->mu = 10000;
	parm_ab->sigma = 5;
	parm_bb->mu = 10;
	parm_bb->sigma = 10;
	parm_cb->mu = 200;
	parm_cb->sigma = 50;
}

int main()
{
	srand(time(NULL));
	testA();
	cout << "start" << endl;
	int statpulls = 31;
	int pulls = 2751;																											// number of pulls, arbitrarily chosen
	agent learner;
	agent* plearner = &learner;
	learner.init();
	arm arm_a;
	arm arm_b;
	arm arm_c;
	arm* parm_a = &arm_a;
	arm* parm_b = &arm_b;
	arm* parm_c = &arm_c;
	arm_a.init();
	arm_b.init();
	arm_c.init();
	vector<double> rewards;
	vector<int> pull_number;
	int arm_a_count = 0;
	int arm_b_count = 0;
	int arm_c_count = 0;
	double action_curve_a = 0;
	double action_curve_b = 0;
	double action_curve_c = 0;
	vector<double> curve_a;
	vector<double> curve_b;
	vector<double> curve_c;
	int q;
	cout << "press 1 to run test b, press any other number to run normal tests" << endl;
	cin >> q;
	if (q == 1) {
		testB(parm_a, parm_b, parm_c);
	}
	cout << "pulling" << endl;
	for (int i = 1; i < statpulls; i++) {
		for (int j = 1; j < pulls; j++) {
			int action_value;
			decide(&learner);																									// learner will choose what to do
			action_value = decide(&learner);																					// places action value into variable to be used later
			double decision_value;																								// 
			decision_value = act1(action_value, parm_a->V_t, parm_b->V_t, parm_c->V_t, parm_a, parm_b, parm_c);					// arm pulled #
			double reward;																										//
			reward = act2(action_value, parm_a->V_t, parm_b->V_t, parm_c->V_t, parm_a, parm_b, parm_c, decision_value);			// reward #
			rewards.push_back(reward);
			double new_expected_value;																							// new expected value
			if (decision_value == 1) {
				new_expected_value = react1(reward, &learner, parm_a, plearner->alpha, parm_a->V_t);
				parm_a->V_t = new_expected_value;
				arm_a_count++;
				curve_a.push_back(action_curve_a++ / j);
			}
			if (decision_value == 2) {
				new_expected_value = react2(reward, &learner, parm_b, plearner->alpha, parm_b->V_t);
				parm_b->V_t = new_expected_value;
				arm_b_count++;
				curve_b.push_back(action_curve_b++ / j);
			}
			if (decision_value == 3) {
				new_expected_value = react3(reward, &learner, parm_c, plearner->alpha, parm_c->V_t);
				parm_c->V_t = new_expected_value;
				arm_c_count++;
				curve_c.push_back(action_curve_c++ / j);
			}
			pull_number.push_back(j);
			double x;
			x = update_alpha(plearner, j);
			plearner->alpha = x;
			/*cout << "alpha after update " << plearner->alpha << endl;*/
			double y;
			y = update_epsilon(plearner, j);
			plearner->epsilon = y;
			/*cout << "epsilon after update " << plearner->epsilon << endl;*/

		}
		cout << "alpha " << plearner->alpha << endl;
		cout << "epsilon " << plearner->epsilon << endl;
		cout << "Arm A pulls: " << arm_a_count << " Arm B pulls: " << arm_b_count << " Arm C pulls: " << arm_c_count << endl;
		cout << "arm a - Mean: " << parm_a->mu << "\t Variance: " << parm_a->sigma << "\t Expected Value: " << parm_a->V_t << endl;
		cout << "arm b - Mean: " << parm_b->mu << "\t Variance: " << parm_b->sigma << "\t Expected Value: " << parm_b->V_t << endl;
		cout << "arm c - Mean: " << parm_c->mu << "\t Variance: " << parm_c->sigma << "\t Expected Value: " << parm_c->V_t << endl;
		cout << "what's happening" << endl;
		ofstream outFile;																										// output file
		outFile.open("Ta_Alton_493_ProjectAlpha.txt");																			// name of output file
		for (int w = 0; w < pulls; w++) {
			outFile << pull_number.at(w) << "\t" << rewards.at(w) << endl;		// outputs pull # and it's corresponding reward to a text file, action curves		
		}
		outFile.close();
	}
	cout << "The program has reached the end and dropped the One Ring into the fires of Mount Doom." << endl;
	return 0;
}