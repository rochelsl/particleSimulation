#pragma once

#include <SFML/Graphics.hpp>
#include <random>
#include <vector>

#include "particle.hpp"

namespace sim {

class Simulation {
public:
    Simulation();

    void step(float dt, sf::Vector2f externalFieldDirection, float externalFieldStrength);

    const std::vector<Particle>& particles() const;

private:
    float temperatureAtTime(float t) const;
    int cellIndex(float x, int n, float L) const;

    sf::Vector2f computeDipoleForceOnP2(sf::Vector2f r, sf::Vector2f mu1, sf::Vector2f mu2, float strength) const;
    sf::Vector2f dipoleField(sf::Vector2f r, sf::Vector2f mu, float strength) const;

    void computeForces();
    void integrateFirstHalf(float dt);
    void integrateSecondHalf(float dt);
    void applyLangevinThermostat(float dt, float temperature);

    std::vector<Particle> particles_;
    std::mt19937 rng_;
    float simulationTime_ = 0.f;

    sf::Vector2f externalFieldDirection_{1.f, 0.2f};
    float externalFieldStrength_ = 100.f;

    float cellSize_ = 1.f;
    int nx_ = 1;
    int ny_ = 1;
    std::vector<std::vector<std::vector<int>>> grid_;
};

}  // namespace sim
