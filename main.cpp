#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include "barrier_option.hpp"


void runTests() { // Test: Up-and-Out + Up-and-In should match vanilla option
    double S0 = 1.0, E = 1.0, r = 0.08, sigma = 0.30, T = 5.0;
    int n_steps = 100;
    double bsPrice = blackScholesCall(S0, E, r, sigma, T);
    double Barrier = 1.5 * S0;

    BarrierOption upAndOut(S0, E, Barrier, r, sigma, T, n_steps, BarrierOption::OptionType::Call, BarrierOption::BarrierStyle::UpAndOut);
    double upAndOutPrice = upAndOut.priceOption(10000);
    BarrierOption upAndIn(S0, E, Barrier, r, sigma, T, n_steps, BarrierOption::OptionType::Call, BarrierOption::BarrierStyle::UpAndIn);
    double upAndInPrice = upAndIn.priceOption(10000);

    double combinedPrice = upAndOutPrice + upAndInPrice;
    assert(std::fabs(combinedPrice - bsPrice) < 0.1);
    std::cout << "Test passed: Up-and-Out + Up-and-In matches vanilla price." << std::endl;
}

void deltaTest() {
    double S0 = 1.0, E = 1.0, r = 0.08, sigma = 0.30, T = 5.0;
    int n_steps = 40, N_sims = 100;
    double h = 0.01; // Increment for finite difference

    BarrierOption upAndOut(S0, E, 2.0 * S0, r, sigma, T, n_steps, BarrierOption::OptionType::Call, BarrierOption::BarrierStyle::UpAndOut);
    std::ofstream deltaupAndOutFile("delta_upAndOut.txt");
    for (double S = 0.0 * S0; S <= 10.0 * S0; S += 0.01) {
        double delta = upAndOut.calculateDelta(S, h, N_sims);
        deltaupAndOutFile << S << " " << delta << std::endl;
    }
    deltaupAndOutFile.close();
    std::cout << "Delta results written to 'delta_upAndOut.txt'." << std::endl;

    BarrierOption upAndIn(S0, E, 9.0 * S0, r, sigma, T, n_steps, BarrierOption::OptionType::Call, BarrierOption::BarrierStyle::UpAndIn);
    std::ofstream deltaupAndInFile("delta_upAndIn.txt");
    for (double S = 0.0 * S0; S <= 10.0 * S0; S += 0.01) {
        double delta = upAndIn.calculateDelta(S, h, N_sims);
        deltaupAndInFile << S << " " << delta << std::endl;
    }
    deltaupAndInFile.close();
    std::cout << "Delta results written to 'delta_upAndIn.txt'." << std::endl;
}

int main() {
    double S0 = 1.0;    // Initial price
    double T = 5.0;      // Time to maturity
    double r = 0.08;     // Risk-free rate
    double sigma = 0.30; // Volatility
    double E = 1.0;      // Strike

    int n_steps = 30;    // number of steps
    int N_sims = 10000;  // number of simulations

    // Part (a): a single path simulation
    {
        BarrierOption demoOption(S0, E, 1.5 * S0, r, sigma, T, n_steps);
        std::mt19937_64 gen(std::random_device{}());
        std::normal_distribution<double> dist(0.0, 1.0);
        auto path = demoOption.simulatePath(gen, dist);

        std::cout << "Sample path:" << std::endl;
        for (auto p : path) {
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }

    // Part (b): Estimate probability of hitting barrier for barriers = 1.5*S0 to 6*S0
    {
        BarrierOption probOption(S0, E, 1.5 * S0, r, sigma, T, n_steps);
        std::cout << "\nProbability of hitting barrier from below:" << std::endl;
        for (double barrierFactor = 1.5; barrierFactor <= 6.0; barrierFactor += 0.5) {
            double barrierLevel = barrierFactor * S0;
            double prob = probOption.estimateBarrierHitProbability(N_sims, barrierLevel);
            std::cout << "Barrier = " << barrierLevel << ", Probability = " << prob << std::endl;
        }
    }

    // Part (c): Using simulation, price the continuous barrier call option
    // For Up-and-Out and Up-and-In, vary N_sims from 10 to 1000 in increments of 10
    {
        std::ofstream out_up_and_out("up_and_out_prices.txt");
        std::ofstream out_up_and_in("up_and_in_prices.txt");

        BarrierOption upAndOutOption(S0, E, 2.0 * S0, r, sigma, T, n_steps, BarrierOption::OptionType::Call, BarrierOption::BarrierStyle::UpAndOut);
        BarrierOption upAndInOption(S0, E, 2.0 * S0, r, sigma, T, n_steps, BarrierOption::OptionType::Call, BarrierOption::BarrierStyle::UpAndIn);

        for (int sims = 10; sims <= 1000; sims += 10) {
            double priceOut = upAndOutOption.priceOption(sims);
            double priceIn = upAndInOption.priceOption(sims);

            out_up_and_out << sims << " " << priceOut << "\n";
            out_up_and_in << sims << " " << priceIn << "\n";
        }

        out_up_and_out.close();
        out_up_and_in.close();

        std::cout << "\nPricing complete. Check 'up_and_out_prices.txt' and 'up_and_in_prices.txt'." << std::endl;
    }
    runTests();
    deltaTest();
    return 0;
}



#endif // MAIN_HPP