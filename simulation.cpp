#include "simulation.hpp"

#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "math_utils.hpp"

namespace sim {

Simulation::Simulation()
    : rng_(std::random_device{}()) {
    cellSize_ = std::max(kLjCutoff, kDipoleCutoff);
    nx_ = std::max(1, static_cast<int>(std::ceil(kWindowWidth / cellSize_)));
    ny_ = std::max(1, static_cast<int>(std::ceil(kWindowHeight / cellSize_)));
    grid_.assign(nx_, std::vector<std::vector<int>>(ny_));

    particles_.reserve(kParticleCount);

    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> angleDist(0.f, 2.f * PI);
    std::normal_distribution<float> normal(0.f, 1.f);

    const int cols = static_cast<int>(std::ceil(std::sqrt(kParticleCount * static_cast<float>(kWindowWidth) / kWindowHeight)));
    const int rows = static_cast<int>(std::ceil(static_cast<float>(kParticleCount) / cols));
    const float dx = static_cast<float>(kWindowWidth) / cols;
    const float dy = static_cast<float>(kWindowHeight) / rows;

    const float initVelSigma = std::sqrt(kInitialTemperature / kMass);
    const float initAngVelSigma = std::sqrt(kInitialTemperature / kMomentOfInertia);

    for (int i = 0; i < kParticleCount; ++i) {
        const int ix = i % cols;
        const int iy = i / cols;

        Particle p;
        p.radius = kParticleRadius;
        p.position = {(ix + 0.5f) * dx + jitter(rng_), (iy + 0.5f) * dy + jitter(rng_)};
        p.velocity = initVelSigma * sf::Vector2f(normal(rng_), normal(rng_));
        p.angle = angleDist(rng_);
        p.angularVelocity = initAngVelSigma * normal(rng_);
        p.magneticMoment = momentFromAngle(p.angle);
        particles_.push_back(p);
    }

    computeForces();
}

void Simulation::step(float dt, sf::Vector2f externalFieldDirection, float externalFieldStrength) {
    externalFieldDirection_ = normalize(externalFieldDirection);
    externalFieldStrength_ = externalFieldStrength;

    const float temperature = temperatureAtTime(simulationTime_);

    integrateFirstHalf(dt);
    computeForces();
    integrateSecondHalf(dt);
    applyLangevinThermostat(dt, temperature);

    simulationTime_ += dt;
}

const std::vector<Particle>& Simulation::particles() const { return particles_; }

float Simulation::temperatureAtTime(float t) const {
    return kFinalTemperature + (kInitialTemperature - kFinalTemperature) * std::exp(-t / kCoolingTime);
}

int Simulation::cellIndex(float x, int n, float L) const {
    x = std::fmod(x, L);
    if (x < 0.f) x += L;

    int c = static_cast<int>(x / cellSize_);
    if (c < 0) c = 0;
    if (c >= n) c = n - 1;
    return c;
}

sf::Vector2f Simulation::computeDipoleForceOnP2(sf::Vector2f r, sf::Vector2f mu1, sf::Vector2f mu2, float strength) const {
    const float r2 = dot(r, r);
    if (r2 < 1e-6f) return {0.f, 0.f};

    const float dist = std::sqrt(r2);
    const sf::Vector2f rHat = r / dist;
    const float effectiveDist = std::max(dist, kDipoleSoftCoreDistance);

    const float mu1DotMu2 = dot(mu1, mu2);
    const float mu1DotR = dot(mu1, rHat);
    const float mu2DotR = dot(mu2, rHat);

    const float prefactor = 3.f * strength / (effectiveDist * effectiveDist * effectiveDist * effectiveDist);
    const sf::Vector2f bracket =
        mu1DotR * mu2 + mu2DotR * mu1 + mu1DotMu2 * rHat - 5.f * mu1DotR * mu2DotR * rHat;

    return prefactor * bracket;
}

sf::Vector2f Simulation::dipoleField(sf::Vector2f r, sf::Vector2f mu, float strength) const {
    const float r2 = dot(r, r);
    if (r2 < 1e-6f) return {0.f, 0.f};

    const float dist = std::sqrt(r2);
    const sf::Vector2f rHat = r / dist;
    const float effectiveDist = std::max(dist, kDipoleSoftCoreDistance);
    const float muDotR = dot(mu, rHat);

    return strength * (3.f * muDotR * rHat - mu) / (effectiveDist * effectiveDist * effectiveDist);
}

void Simulation::computeForces() {
    for (auto& col : grid_) {
        for (auto& cell : col) cell.clear();
    }

    const sf::Vector2f Bext = externalFieldStrength_ * externalFieldDirection_;

    for (auto& p : particles_) {
        applyPBC(p.position, static_cast<float>(kWindowWidth), static_cast<float>(kWindowHeight));
        p.acceleration = {0.f, 0.f};
        p.angularAcceleration = cross2D(p.magneticMoment, Bext) / kMomentOfInertia;
    }

    for (int i = 0; i < static_cast<int>(particles_.size()); ++i) {
        const int cx = cellIndex(particles_[i].position.x, nx_, static_cast<float>(kWindowWidth));
        const int cy = cellIndex(particles_[i].position.y, ny_, static_cast<float>(kWindowHeight));
        grid_[cx][cy].push_back(i);
    }

    const float ljCutoff2 = kLjCutoff * kLjCutoff;
    const float dipoleCutoff2 = kDipoleCutoff * kDipoleCutoff;
    const float sr2c = (kSigma * kSigma) / ljCutoff2;
    const float sr6c = sr2c * sr2c * sr2c;
    const float sr12c = sr6c * sr6c;
    const float ljScalarAtCutoff = 24.f * kEpsilon * (2.f * sr12c - sr6c) / ljCutoff2;

    for (int cx = 0; cx < nx_; ++cx) {
        for (int cy = 0; cy < ny_; ++cy) {
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    const int cx2 = (cx + dx + nx_) % nx_;
                    const int cy2 = (cy + dy + ny_) % ny_;

                    for (int i : grid_[cx][cy]) {
                        for (int j : grid_[cx2][cy2]) {
                            if (i >= j) continue;

                            Particle& p1 = particles_[i];
                            Particle& p2 = particles_[j];

                            sf::Vector2f r = minimumImage(p2.position - p1.position, static_cast<float>(kWindowWidth), static_cast<float>(kWindowHeight));
                            const float r2 = dot(r, r);
                            if (r2 < 1e-6f) continue;

                            sf::Vector2f totalForceOnP2{0.f, 0.f};

                            if (r2 < ljCutoff2) {
                                const float invR2 = 1.f / r2;
                                const float sr2 = (kSigma * kSigma) * invR2;
                                const float sr6 = sr2 * sr2 * sr2;
                                const float sr12 = sr6 * sr6;

                                const float ljScalar = 24.f * kEpsilon * invR2 * (2.f * sr12 - sr6) - ljScalarAtCutoff;
                                totalForceOnP2 += r * ljScalar;
                            }

                            if (r2 < dipoleCutoff2) {
                                totalForceOnP2 += computeDipoleForceOnP2(r, p1.magneticMoment, p2.magneticMoment, kDipoleStrength);

                                const sf::Vector2f BOn1 = dipoleField(neg(r), p2.magneticMoment, kDipoleStrength);
                                const sf::Vector2f BOn2 = dipoleField(r, p1.magneticMoment, kDipoleStrength);

                                p1.angularAcceleration += cross2D(p1.magneticMoment, BOn1) / kMomentOfInertia;
                                p2.angularAcceleration += cross2D(p2.magneticMoment, BOn2) / kMomentOfInertia;
                            }

                            p1.acceleration -= totalForceOnP2 / kMass;
                            p2.acceleration += totalForceOnP2 / kMass;
                        }
                    }
                }
            }
        }
    }
}

void Simulation::integrateFirstHalf(float dt) {
    for (auto& p : particles_) {
        p.position += p.velocity * dt + p.acceleration * (0.5f * dt * dt);
        p.velocity += p.acceleration * (0.5f * dt);

        p.angle += p.angularVelocity * dt + 0.5f * p.angularAcceleration * dt * dt;
        p.angularVelocity += 0.5f * p.angularAcceleration * dt;

        p.magneticMoment = momentFromAngle(p.angle);
        applyPBC(p.position, static_cast<float>(kWindowWidth), static_cast<float>(kWindowHeight));
    }
}

void Simulation::integrateSecondHalf(float dt) {
    for (auto& p : particles_) {
        p.velocity += p.acceleration * (0.5f * dt);
        p.angularVelocity += 0.5f * p.angularAcceleration * dt;
        p.magneticMoment = momentFromAngle(p.angle);
    }
}

void Simulation::applyLangevinThermostat(float dt, float temperature) {
    std::normal_distribution<float> normal(0.f, 1.f);

    const float cT = std::exp(-kGammaT * dt);
    const float cR = std::exp(-kGammaR * dt);

    const float sigmaV = std::sqrt((temperature / kMass) * (1.f - cT * cT));
    const float sigmaW = std::sqrt((temperature / kMomentOfInertia) * (1.f - cR * cR));

    for (auto& p : particles_) {
        p.velocity = cT * p.velocity + sigmaV * sf::Vector2f(normal(rng_), normal(rng_));
        p.angularVelocity = cR * p.angularVelocity + sigmaW * normal(rng_);
    }
}

}  // namespace sim
