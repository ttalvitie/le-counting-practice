#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <cmath>
#include <vector>
#include <cstdio>

#include "tools.hpp"


// the Stopping Rule Algorithm from Dagum et al., Optimal Monte Carlo Estimation,
// produces a Monte Carlo estimate of the mean of a Bernoulli random variable with given
// precision guarantees after being given a sufficient number of samples
struct StoppingRuleEstimator
{
	int n_samples;  // number of samples given
	double sum;     // cumulative sum of smaples
	double Y1;      // required number of samples to give an estimate with precision guarantees
	
	// sets the desired precision guarantees, an (ε,δ)-estimate
	// this can be changed at any time (the samples given so far are not reset)
	void set_precision(double epsilon, double delta)
	{
		double lambda = exp(1) - 2;
		double Y = 4 * lambda * log(2 / delta) / pow(epsilon, 2);
		Y1 = 1 + (1 + epsilon) * Y;
	}
	
	// resets the sampler to the initial state where no samples have been seen
	void reset()
	{
		n_samples = 0.0;
		sum = 0.0;
	}
	
	// initializes the estimator without precision guarantees (must be set before asking for estimates)
	StoppingRuleEstimator()
	{
		reset();
	}
	
	// initializes the estimator with precision guarantees
	StoppingRuleEstimator(double epsilon, double delta)
	{
		reset();
		set_precision(epsilon, delta);
	}
	
	// feeds a new sample to the estimator
	void feed(double x)
	{
		++n_samples;
		sum += x;
	}
	
	// whether enough samples to give an estimate with the specified precision guarantees
	bool is_sufficient()
	{
		return sum >= Y1;
	}
	
	// Monte Carlo estimate based on samples so far, no precision guarantees
	double preliminary_estimate()
	{
		return sum / n_samples;
	}
	
	// Monte Carlo estimate with the specified precision guarantees once is_sufficient() is true
	double estimate()
	{
		return Y1 / n_samples;
	}
};

#endif
