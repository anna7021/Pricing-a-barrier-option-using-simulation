#include "barrier_option.hpp"
#include <cmath>
#include <algorithm>

BarrierOption::BarrierOption(double S0, double E, double B, double r, double sigma, double T, int n_steps,
    OptionType optType, BarrierStyle style)
    : m_S0(S0), m_E(E), m_B(B), m_r(r), m_sigma(sigma), m_T(T), m_n_steps(n_steps),
    m_optType(optType), m_style(style)
{
}

std::vector<double> BarrierOption::simulatePath(std::mt19937_64& gen, std::normal_distribution<double>& dist, bool invert) const {
    std::vector<double> path(m_n_steps + 1, 0.0);
    path[0] = m_S0;
    double dt = m_T / static_cast<double>(m_n_steps);
    double drift = (m_r - 0.5 * m_sigma * m_sigma) * dt;
    double vol = m_sigma * std::sqrt(dt);

    for (int i = 1; i <= m_n_steps; ++i) {
        double Z = dist(gen);
        if (invert) Z = -Z; // Use antithetic variate
        path[i] = path[i - 1] * std::exp(drift + vol * Z);
    }
    return path;
}

bool BarrierOption::barrierHit(const std::vector<double>& path) const {
    // for "up" type barriers, we consider if max element of path >= B
    for (auto& price : path) {
        if (price >= m_B) return true;
    }
    return false;
}

double BarrierOption::payoff(const std::vector<double>& path) const {
    double S_final = path.back();
    double intrinsic = 0.0;
    if (m_optType == OptionType::Call)
        intrinsic = std::max(S_final - m_E, 0.0);
    else
        intrinsic = std::max(m_E - S_final, 0.0);

    bool hit = barrierHit(path);

    if (m_style == BarrierStyle::UpAndOut) { // For Up-And-Out: payoff is intrinsic if barrier not hit, else 0
        if (hit) return 0.0;
        return intrinsic;
    }
    
    if (m_style == BarrierStyle::UpAndIn) { // For Up-And-In: payoff is intrinsic if barrier hit, else 0
        if (hit) return intrinsic;
        return 0.0;
    }

    return 0.0; // default (should not happen)
}

void BarrierOption::setBarrier(double B) {
    m_B = B;
}

void BarrierOption::setBarrierStyle(BarrierStyle style) {
    m_style = style;
}

void BarrierOption::setOptionParameters(double S0, double E, double T, double r, double sigma, int n_steps) {
    m_S0 = S0; m_E = E; m_T = T; m_r = r; m_sigma = sigma; m_n_steps = n_steps;
}

double BarrierOption::estimateBarrierHitProbability(int N_sims, double barrierToTest) {
    double oldBarrier = m_B;
    m_B = barrierToTest;

    std::mt19937_64 gen(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);

    int hitCount = 0;
    for (int i = 0; i < N_sims; ++i) {
        auto path = simulatePath(gen, dist);
        if (barrierHit(path)) {
            hitCount++;
        }
    }
    m_B = oldBarrier;
    return static_cast<double>(hitCount) / N_sims;
}

double BarrierOption::priceOption(int N_sims) {
    std::mt19937_64 gen(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);

    double sumPayoffs = 0.0;
    for (int i = 0; i < N_sims; ++i) {
        auto path1 = simulatePath(gen, dist);
        auto path2 = simulatePath(gen, dist, true); // Use antithetic variate (negative Z)
        sumPayoffs += payoff(path1);
        sumPayoffs += payoff(path2);
    }

    // Discount payoff and average over 2 * N_sims
    double discountFactor = std::exp(-m_r * m_T);
    double price = (sumPayoffs / (2 * N_sims)) * discountFactor;
    return price;
}

double BarrierOption::calculateDelta(double S0, double h, int N_sims) {
    double deltaSum = 0.0;
    int n = 5; // Number of steps for smoothing

    for (int i = 1; i <= n; ++i) {
        double step = i * h;
        this->setOptionParameters(S0 + step, m_E, m_T, m_r, m_sigma, m_n_steps);
        double pricePlusH = this->priceOption(N_sims);

        this->setOptionParameters(S0 - step, m_E, m_T, m_r, m_sigma, m_n_steps);
        double priceMinusH = this->priceOption(N_sims);

        deltaSum += (pricePlusH - priceMinusH) / (2 * step);
    }

    this->setOptionParameters(S0, m_E, m_T, m_r, m_sigma, m_n_steps); // Restore original parameters

    return deltaSum / n; 
}

double blackScholesCall(double S0, double E, double r, double sigma, double T) {
    double d1 = (std::log(S0 / E) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    double N_d1 = 0.5 * std::erfc(-d1 / std::sqrt(2.0)); // complementary error function
    double N_d2 = 0.5 * std::erfc(-d2 / std::sqrt(2.0)); 
    return S0 * N_d1 - E * std::exp(-r * T) * N_d2;
}
