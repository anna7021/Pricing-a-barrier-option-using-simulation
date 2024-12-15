#ifndef BARRIER_OPTION_HPP
#define BARRIER_OPTION_HPP

#include <vector>
#include <random>
#include <cmath>

class BarrierOption {
public:
	enum class OptionType {Call, Put};
	enum class BarrierStyle{UpAndOut, UpAndIn};
	BarrierOption(double S0, double E, double B, double r, double sigma, double T, int n_steps,
				OptionType optType = OptionType::Call, BarrierStyle style = BarrierStyle::UpAndOut);

	// Simulate a single asset price path under GBM using n_steps
	// Returns a vector of prices at each time step (including initial)
	std::vector<double> simulatePath(std::mt19937_64& gen, std::normal_distribution<double>& dist, bool invert = false) const;
	
	bool barrierHit(const std::vector<double>& path) const;   // Check if the barrier was hit at any point in a given path
	double payoff(const std::vector<double>& path) const;	// Compute the payoff of the option given a path
	double estimateBarrierHitProbability(int N_sims, double barrierToTest);	 // Estimate the probability of hitting the barrier from below by simulation
	double priceOption(int N_sims);	 // Price the barrier option using Monte Carlo simulation
	
	void setBarrier(double B);	// Setters to easily change parameters
	void setBarrierStyle(BarrierStyle style);	
	void setOptionParameters(double S0, double E, double T, double r, double sigma, int n_steps);
	double calculateDelta(double S0, double h, int N_sims);

private:
	double m_S0;       // Initial price
	double m_E;        // Strike price
	double m_B;        // Barrier
	double m_r;        // Risk-free rate
	double m_sigma;    // Volatility
	double m_T;        // Time to maturity
	int m_n_steps;     // Number of steps
	OptionType m_optType;
	BarrierStyle m_style;

};

double blackScholesCall(double S0, double E, double r, double sigma, double T);

#endif // BARRIER_OPTION_HPP