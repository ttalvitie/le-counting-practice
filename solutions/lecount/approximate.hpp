#ifndef APPROXIMATE_H
#define APPROXIMATE_H

#include <utility>
#include <cmath>

#include "relaxation.hpp"
#include "montecarlo.hpp"


struct ApproximationOptions
{
	RelaxationOptions relaxation;
	
	double epsilon;
	double delta;
	
	int partition_initial_size;
	int partition_size_increment;
	double tentative_sampling_time;
	double dp_time_scaling;
	
	bool verbose;
	
	ApproximationOptions()
	{
		epsilon = 0.5;
		delta = 0.25;
		
		partition_initial_size = 20;
		partition_size_increment = 5;
		tentative_sampling_time = 0.1;
		dp_time_scaling = 10;
		
		verbose = true;
	}
};



template <typename N>
struct MonteCarloEstimator
{
	ApproximationOptions &options;
	const digraph &poset;
	RelaxationSampler<N> *relaxation;
	StoppingRuleEstimator sre;
	double sampler_creation_time;
	
	int search_relaxation()
	{
		RelaxationOptions &opt = options.relaxation;
		
		if (options.verbose) eprintf("============================================================= Relaxation search, max set size: %i\n", opt.partition_max_size);
		
		// find a new relaxation and try to build the sampler
		double start_time = seconds();
		RelaxationSampler<N> *new_relaxation = new RelaxationSampler<N>(opt, poset);
		bool success = new_relaxation->init_sampler();
		sampler_creation_time = seconds() - start_time;
		
		// if ran out of memory
		if (!success) {
			if (options.verbose) eprintf("Failed to compute the relaxation.\n")
			return 0;
		}
		
		if (options.verbose) eprintf("\nTime: %.2f  Extensions: %Le\n", sampler_creation_time, (long double)new_relaxation->n_orders);
		
		// if there's a previous relaxation and it's better than the new one
		if (relaxation != NULL && relaxation->n_orders <= new_relaxation->n_orders) {
			if (options.verbose) eprintf("Keeping previous relaxation.\n");
			delete new_relaxation;
			return 2;
		}
		
		// otherwise replace the old relaxation (if any) with the new one
		delete relaxation;
		relaxation = new_relaxation;
		return 1;
	}
	
	void choose_relaxation()
	{
		if (options.relaxation.relaxation == RELAXATION_MANUAL) {
			while (search_relaxation() == 0);
			return;
		}
		
		options.relaxation.verbose = options.verbose;
		options.relaxation.partition_max_size = options.partition_initial_size;
		search_relaxation();
		
		while (!evaluate_relaxation()) {
			// if already reached maximum size, accept this anyway
			if (options.relaxation.partition_max_size == poset.n) break;
			
			options.relaxation.partition_max_size += options.partition_size_increment;
			if (options.relaxation.partition_max_size >= poset.n) options.relaxation.partition_max_size = poset.n;
			
			int ret = search_relaxation();
			if (ret == 0) {
				if (options.verbose) eprintf("Settling for the previous relaxation.\n")
				break;
			}
		}
	}
	
	MonteCarloEstimator(ApproximationOptions &options, const digraph &poset) : options(options), poset(poset), relaxation(NULL)
	{
		sre.set_precision(options.epsilon, options.delta);
		if (poset.n != 0) choose_relaxation();
	}
	
	~MonteCarloEstimator()
	{
		delete relaxation;
	}
	
	bool evaluate_relaxation()
	{
		int n_accepts_required = ceil(sre.Y1);
		
		double test_proportion = options.tentative_sampling_time;
		int n_test_accepts_required = n_accepts_required * test_proportion;
		
		double max_test_time = sampler_creation_time * test_proportion;
		
		if (options.verbose) eprintf("Test sampling until %i accepts or %.2f seconds used.\n", n_test_accepts_required, max_test_time);
		
		double start_time = seconds();
		
		int n_samples = 0;
		int n_accepts = 0;
		
		while (seconds() - start_time < max_test_time) {
			if (n_accepts == n_test_accepts_required) break;
			++n_samples;
			if (relaxation->draw_sample()) ++n_accepts;
		}
		
		double test_time = seconds() - start_time;
		
		N ar_estimate = (double)n_accepts / n_samples;
		double estimated_accept_speed = (double)n_accepts / test_time;
		
		if (options.verbose) eprintf("\rAccepted: %i / %i (%.2Le)   Time: %.2f   Accept speed: %.2f\n", n_accepts, n_samples, (long double)ar_estimate, test_time, estimated_accept_speed);
		
		double estimated_time_if_accept = n_accepts_required / estimated_accept_speed;
		double estimated_time_if_reject = sampler_creation_time * options.dp_time_scaling;
		
		if (options.verbose) eprintf("Estimated time on accept/reject: %.2f / %.2f\n", estimated_time_if_accept, estimated_time_if_reject);
		
		return estimated_time_if_accept < estimated_time_if_reject;
	}
	
	N approximate()
	{
		if (poset.n == 0) return 1;
		
		if (options.verbose) {
			eprintf("============================================================= Sampling phase\n");
			relaxation->display_relaxation();
			eprintf("ε = %.4f   δ = %.4f   Y1 = %f\n", options.epsilon, options.delta, sre.Y1);
		}
		
		int n_accepts_required = ceil(sre.Y1);
		
		double start_time = seconds();
		
		// sample until enough accepted linear extensions
		while (!sre.is_sufficient()) {
			int result = relaxation->draw_sample() ? 1 : 0;
			sre.feed(result);
			
			if ((sre.n_samples % 100000 == 0) || sre.is_sufficient()) {
				N ar_estimate = sre.sum / sre.n_samples;
				double time_passed = seconds() - start_time;
				int sample_speed = sre.n_samples / time_passed;
				N accept_speed_estimate = sre.sum / time_passed;
				if (options.verbose) eprintf("\rSamples: %i (%i/s)   Accepted: %i/%i   Rate: %.2Le   Estimate: %.2Le   ETL: %is  ",
					sre.n_samples,
					sample_speed,
					(int)sre.sum,
					n_accepts_required,
					(long double)ar_estimate,
					(long double)ar_estimate * relaxation->n_orders,
					accept_speed_estimate > 0 ? (int)((n_accepts_required - sre.sum) / accept_speed_estimate) : -1);
				eflush();
			}
		}
		
		if (options.verbose) eprintf("\n");
		
		// (ε,δ)-estimate of the acceptance rate
		N acceptance_rate_estimate = sre.estimate();
		
		// (ε,δ)-estimate of the number of linear extensions
		N estimate = acceptance_rate_estimate * relaxation->n_orders;
		
		if (options.verbose) eprintf("Estimate: %Lf\n", (long double)estimate);
		
		return estimate;
	}
};

#endif
