#include <SFML/Graphics.hpp> 
#include <iostream>
#include <vector>
#include <random>
#include <cmath>

constexpr float PI = 3.14159265359f;

//Window size
const int width = 1920;
const int height = 1080;

//Initialize parameters for LJ-potential, cutoff and grid-calculations
const float epsilon = 1.0f;
const float sigma = 10.0f;
const float rc = 1.0f * sigma;
const float rc2 = rc * rc;

//for quantitative real NP behavior
constexpr double metersPerPixel = 1e-9; // 1 pixel = 1 nm 
float radiusPixels = 5.0f; //radius of the nanoparticles

constexpr double kB   = 1.380649e-23;        // J/K
constexpr double mu0  = 4.0 * PI * 1e-7;     // N/A^2
constexpr double muB  = 9.2740100783e-24;    // A m^2

constexpr double T = 300.0;                  // K
constexpr double eta = 1.0e-3;               // Pa s, water

constexpr double particleDiameter = 10e-9;   // 10 nm
constexpr double radius = particleDiameter / 2.0;

constexpr double magneticMoment = 1e5 * muB; // A m^2
//for a sphere in water
constexpr double gammaTrans = 6.0 * PI * eta * radius;
//rotational drag
constexpr double gammaRot = 8.0 * PI * eta * radius * radius * radius;

//dipole strength
const float dipoleStrength = 100.0f;
std::uniform_real_distribution<float> angleDist(0.f, 2.f * PI);
//to be deactivated when Maxwe--Boltzmann temperature is included
//const float rotationalDamping = 0.999; // To decrease probability of numerical explosion
//Background field
sf::Vector2f normalize(sf::Vector2f v) {
    float len = std::sqrt(v.x * v.x + v.y * v.y);
    if (len < 1e-8f)
        return {0.f, 0.f};
    return v / len;
}
const sf::Vector2f externalFieldDirection = normalize({1.0f, 0.0f});
const float externalFieldStrength = 100.0f;

//for temperature and Langevin behavior
const float temperature = 0.2f;
const float gamma = 0.5f;
const float gammaRot = 0.5f;
const float mass = 1.0f;

//For grid list, performance boost
float cellSize = rc;
int nx = static_cast<int>(width / cellSize);
int ny = static_cast<int>(height / cellSize);
std::vector<std::vector<std::vector<int>>> grid(nx, std::vector<std::vector<int>>(ny));

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;
    //magnetic moments
    sf::Vector2f magneticMoment;
    float angle;
    float angularVelocity;
    float angularAcceleration;

    float radius;
    //temperature and Langevin behavior
    float mass;
    float momentOfInertia;

    //for the Brownian translational step of NPs moving in water, the force acting on the particle needs to be stored
    float torque;
    sf::Vector2f force;
};

//random magnetization vector orientation initialization
sf::Vector2f randomUnitVector(std::mt19937& rng) {
    std::uniform_real_distribution<float> angleDist(0.f, 2.f * PI);

    float theta = angleDist(rng);

    return {
        std::cos(theta),
        std::sin(theta)
    };
}
//For the dot product in the dipolar interaction
float dot(sf::Vector2f a, sf::Vector2f b) {
    return a.x * b.x + a.y * b.y;
}

//For the vector product
float cross2D(sf::Vector2f a, sf::Vector2f b) {
    return a.x * b.y - a.y * b.x;
}

sf::Vector2f momentFromAngle(float angle) {
    return {std::cos(angle), std::sin(angle)};
}

sf::Vector2f computeDipoleForceOnP2(
    sf::Vector2f r,
    sf::Vector2f mu1,
    sf::Vector2f mu2,
    float dipoleStrength
) {
    float r2 = dot(r, r);

    if (r2 < 1e-6f)
        return {0.f, 0.f};

    float dist = std::sqrt(r2);
    sf::Vector2f rHat = r / dist;

    float mu1_dot_mu2 = dot(mu1, mu2);
    float mu1_dot_r = dot(mu1, rHat);
    float mu2_dot_r = dot(mu2, rHat);

    float prefactor = 3.f * dipoleStrength / (dist * dist * dist * dist);

    sf::Vector2f bracket =
        mu1_dot_r * mu2 +
        mu2_dot_r * mu1 +
        mu1_dot_mu2 * rHat -
        5.f * mu1_dot_r * mu2_dot_r * rHat;

    return prefactor * bracket; // force on particle 2
}

sf::Vector2f dipoleField(
    sf::Vector2f r,
    sf::Vector2f mu,
    float dipoleStrength
) {
    float r2 = dot(r, r);

    if (r2 < 1e-6f)
        return {0.f, 0.f};

    float dist = std::sqrt(r2);
    sf::Vector2f rHat = r / dist;

    float muDotR = dot(mu, rHat);

    return dipoleStrength *
           (3.f * muDotR * rHat - mu) /
           (dist * dist * dist);
}

// OBSELETE NOW THAT THE PBC IS INCLUDED
// void handleWallCollision(Particle& p, const sf::Vector2u& size) {
//     if (p.position.x - p.radius < 0.f) {
//         p.position.x = p.radius;
//         p.velocity.x *= -1.f;
//     } else if (p.position.x + p.radius > size.x) {
//         p.position.x = size.x - p.radius;
//         p.velocity.x *= -1.f;
//     }

//     if (p.position.y - p.radius < 0.f) {
//         p.position.y = p.radius;
//         p.velocity.y *= -1.f;
//     } else if (p.position.y + p.radius > size.y) {
//         p.position.y = size.y - p.radius;
//         p.velocity.y *= -1.f;
//     }
// }

void applyPBC(sf::Vector2f& pos, float w, float h) {
    pos.x = std::fmod(pos.x, w);
    pos.y = std::fmod(pos.y, h);

    if (pos.x < 0.f) pos.x += w;
    if (pos.y < 0.f) pos.y += h;
}

sf::Vector2f minimumImage(sf::Vector2f r, float w, float h) {
    if (r.x >  w / 2.f) r.x -= w;
    if (r.x < -w / 2.f) r.x += w;

    if (r.y >  h / 2.f) r.y -= h;
    if (r.y < -h / 2.f) r.y += h;

    return r;
}

int cellIndex(float x, int n, float L) {
    x = std::fmod(x, L);
    if (x < 0.f) x += L;

    int c = static_cast<int>(x / cellSize);
    if (c < 0) c = 0;
    if (c >= n) c = n - 1;
    return c;
}

void computeForces(std::vector<Particle>& particles) {
    // Clear grid at the beginning, not at the end.
    for (auto& col : grid)
        for (auto& cell : col)
            cell.clear();

    // Reset accelerations and wrap positions before cell assignment.
    for (auto& p : particles) {
        p.acceleration = {0.f, 0.f};
        p.angularAcceleration = 0.f; //magnetic moment angular acceleration
        applyPBC(p.position, width, height);
    }

    // Assign particles to cells.
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        int cx = cellIndex(particles[i].position.x, nx, width);
        int cy = cellIndex(particles[i].position.y, ny, height);
        grid[cx][cy].push_back(i);
    }

    // Force-shifted LJ cutoff value. This is the scalar multiplying r.
    const float sr2c = (sigma * sigma) / rc2;
    const float sr6c = sr2c * sr2c * sr2c;
    const float sr12c = sr6c * sr6c;
    const float f_shift = 24.f * epsilon * (2.f * sr12c - sr6c) / rc2;

    for (int cx = 0; cx < nx; ++cx) {
        for (int cy = 0; cy < ny; ++cy) {
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    int cx2 = (cx + dx + nx) % nx;
                    int cy2 = (cy + dy + ny) % ny;

                    for (int i : grid[cx][cy]) {
                        for (int j : grid[cx2][cy2]) {
                            if (i >= j) continue;

                            Particle& p1 = particles[i];
                            Particle& p2 = particles[j];

                            sf::Vector2f r = p2.position - p1.position;
                            r = minimumImage(r, width, height);

                            //For magnetization vector
                            sf::Vector2f B_on_1 = dipoleField(-r, p2.magneticMoment, dipoleStrength);
                            sf::Vector2f B_on_2 = dipoleField( r, p1.magneticMoment, dipoleStrength);
                            float torque1 = cross2D(p1.magneticMoment, B_on_1);
                            float torque2 = cross2D(p2.magneticMoment, B_on_2);
                            //Background magnetic field
                            sf::Vector2f B_ext = externalFieldStrength * externalFieldDirection;
                            for (auto& p : particles) {
                                float externalTorque = cross2D(p.magneticMoment, B_ext);
                                //magnetic field only
                                //p.angularAcceleration += externalTorque;
                                //magnetic field and temperature
                                p.angularAcceleration += externalTorque / p.momentOfInertia;
                            }
                            //p1.angularAcceleration += torque1; //This assumed moment of inertia = 1: angularAcceleration += torque / 1
                            //p2.angularAcceleration += torque2;
                            p1.angularAcceleration += torque1 / p1.momentOfInertia;
                            p2.angularAcceleration += torque2 / p2.momentOfInertia;

                            float r2 = r.x * r.x + r.y * r.y;
                            if (r2 > rc2 || r2 < 1e-6f) continue;

                            // Lennard-Jones force
                            float inv_r2 = 1.f / r2;
                            float sr2 = (sigma * sigma) * inv_r2;
                            float sr6 = sr2 * sr2 * sr2;
                            float sr12 = sr6 * sr6;

                            float ljScalar = 24.f * epsilon * inv_r2 * (2.f * sr12 - sr6);

                            sf::Vector2f ljForce = r * ljScalar;

                            // Dipolar force
                            sf::Vector2f dipoleForce = computeDipoleForceOnP2(
                                r,
                                p1.magneticMoment,
                                p2.magneticMoment,
                                dipoleStrength
                            );

                            sf::Vector2f totalForce = ljForce + dipoleForce;

                            // Newton's third law
                            //p1.acceleration -= totalForce;
                            //p2.acceleration += totalForce;
                            p1.acceleration -= totalForce / p1.mass;
                            p2.acceleration += totalForce / p2.mass;
                        }
                    }
                }
            }
        }
    }
}

//Without magnetization vector rotational integration
// void integrate(std::vector<Particle>& particles, float dt) {
//     for (auto& p : particles) {
//         p.position += p.velocity * dt + 0.5f * p.acceleration * dt * dt;
//         applyPBC(p.position, width, height);
//         p.velocity += 0.5f * p.acceleration * dt;
//     }
// }

void integrate(std::vector<Particle>& particles, float dt) {
    for (auto& p : particles) {
        p.position += p.velocity * dt + 0.5f * p.acceleration * dt * dt;
        p.velocity += 0.5f * p.acceleration * dt;

        p.angle += p.angularVelocity * dt
                 + 0.5f * p.angularAcceleration * dt * dt;

        p.angularVelocity += 0.5f * p.angularAcceleration * dt;

        p.magneticMoment = momentFromAngle(p.angle);
    }
}

void applyLangevinThermostat(
    std::vector<Particle>& particles,
    float dt,
    float temperature,
    float gamma,
    std::mt19937& rng
) {
    std::normal_distribution<float> normal(0.f, 1.f);

    for (auto& p : particles) {
        float noiseScale = std::sqrt(2.f * gamma * temperature * dt / p.mass);

        p.velocity.x += -gamma * p.velocity.x * dt + noiseScale * normal(rng);
        p.velocity.y += -gamma * p.velocity.y * dt + noiseScale * normal(rng);
    }
}

void applyRotationalLangevinThermostat(
    std::vector<Particle>& particles,
    float dt,
    float temperature,
    float gammaRot,
    std::mt19937& rng
) {
    std::normal_distribution<float> normal(0.f, 1.f);

    for (auto& p : particles) {
        float noiseScale = std::sqrt(
            2.f * gammaRot * temperature * dt / p.momentOfInertia
        );

        p.angularVelocity += -gammaRot * p.angularVelocity * dt
                           + noiseScale * normal(rng);
    }
}

void brownianStep(
    std::vector<Particle>& particles,
    double dt,
    std::mt19937& rng
) {
    std::normal_distribution<double> normal(0.0, 1.0);

    const double Dt = kB * T / gammaTrans;
    const double noiseScale = std::sqrt(2.0 * Dt * dt);

    for (auto& p : particles) {
        p.position.x += static_cast<float>((p.force.x / gammaTrans) * dt
                      + noiseScale * normal(rng));

        p.position.y += static_cast<float>((p.force.y / gammaTrans) * dt
                      + noiseScale * normal(rng));
    }
}

void brownianRotationStep(
    std::vector<Particle>& particles,
    double dt,
    std::mt19937& rng
) {
    std::normal_distribution<double> normal(0.0, 1.0);

    const double Dr = kB * T / gammaRot;
    const double noiseScale = std::sqrt(2.0 * Dr * dt);

    for (auto& p : particles) {
        p.angle += static_cast<float>((p.torque / gammaRot) * dt
                 + noiseScale * normal(rng));

        p.magneticMoment = momentFromAngle(p.angle);
    }
}

std::vector<float> computeRadialDistribution(
    const std::vector<Particle>& particles,
    float width,
    float height,
    int bins,
    float rMax
) {
    std::vector<float> hist(bins, 0.f);

    const int N = static_cast<int>(particles.size());
    const float dr = rMax / bins;
    const float area = width * height;
    const float rho = N / area;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            sf::Vector2f r = particles[j].position - particles[i].position;
            r = minimumImage(r, width, height);

            float dist = std::sqrt(r.x * r.x + r.y * r.y);

            if (dist < rMax) {
                int bin = static_cast<int>(dist / dr);
                hist[bin] += 2.f; // count both i-j and j-i
            }
        }
    }

    std::vector<float> g(bins, 0.f);

    for (int k = 0; k < bins; ++k) {
        float rInner = k * dr;
        float rOuter = rInner + dr;

        float shellArea = static_cast<float>(M_PI) *
                          (rOuter * rOuter - rInner * rInner);

        float idealCount = N * rho * shellArea;

        if (idealCount > 0.f)
            g[k] = hist[k] / idealCount;
    }

    return g;
}

int main() { 
    auto window = sf::RenderWindow{{width, height}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    std::vector<Particle> particles;

    const int N = 1000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    //for random distribution of velocities
    //std::uniform_real_distribution<float> vel(-10.f, 10.f);
    //for Maxwell-Boltzman-like distribution of velocities
    std::normal_distribution<float> velDist(0.f, std::sqrt(temperature / mass));


    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = radiusPixels;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        //for random distribution of velocities
        //p.velocity = {vel(rng), vel(rng)};
        //p.acceleration = {0.f, 0.f};
        //for temperature and Langevin
        p.velocity = {velDist(rng), velDist(rng)};
        //magnetic moment
        p.magneticMoment = randomUnitVector(rng); //magnetization vector
        p.angle = angleDist(rng);
        //for random angular velocities
        //p.angularVelocity = 10.f;
        //for Maxwell-Boltzmann-like velocities
        std::normal_distribution<float> omegaDist(
            0.f,
            std::sqrt(temperature / p.momentOfInertia)
        );
        p.angularVelocity = omegaDist(rng);
        p.angularAcceleration = 0.f;
        p.magneticMoment = momentFromAngle(p.angle);
        //temperature and Langevin behavior
        p.momentOfInertia = 0.5f * p.mass * p.radius * p.radius;
        particles.push_back(p);
    }

    // Initial acceleration for velocity Verlet.
    computeForces(particles);

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    //sf::Clock clock;
    //float dt = clock.restart().asSeconds();
    const float dt = 0.001f;

    //For radial distribution function
    // int frameCounter = 0;
    // std::vector<float> gr;
    // const int grBins = 200;
    // const float rMax = std::min(width, height) * 0.5f;

    // std::vector<float> grAverage(grBins, 0.f);
    // int grSamples = 0;

    while (window.isOpen()) { 

        for (auto event = sf::Event{}; window.pollEvent(event);) { 
            if (event.type == sf::Event::Closed) { 
                window.close(); 
            }
			else if (event.type == sf::Event::KeyPressed) {
				if (event.key.code == sf::Keyboard::Escape) {
					window.close();
				}
			}
        }

        // 1. integrate positions + half velocity
        integrate(particles, dt);
        // 2. recompute forces
        computeForces(particles);
        // 3. finish velocity update
        for (auto& p : particles) {
            p.velocity += 0.5f * p.acceleration * dt;
            p.angularVelocity += 0.5f * p.angularAcceleration * dt;

            p.magneticMoment = momentFromAngle(p.angle);

            applyPBC(p.position, width, height);
        }

        applyLangevinThermostat(particles, dt, temperature, gamma, rng);
        applyRotationalLangevinThermostat(particles, dt, temperature, gammaRot, rng);
        // 4. walls
        // for (auto& p : particles)
        //     handleWallCollision(p, window.getSize());

        window.clear();
        //Collision with the wall check
        for (const auto& p : particles) {
            shape.setRadius(p.radius);
            shape.setOrigin({p.radius, p.radius});
            shape.setPosition(p.position);
            window.draw(shape);

            // Draw the magnetization axis
            sf::Vertex line[] = {
            sf::Vertex(p.position, sf::Color::Red),
            sf::Vertex(p.position + p.magneticMoment * p.radius * 2.0f, sf::Color::Red)
        };
        window.draw(line, 2, sf::Lines);
        }

        sf::Vertex fieldLine[] = {
            sf::Vertex(sf::Vector2f(50.f, 50.f), sf::Color::Blue),
            sf::Vertex(sf::Vector2f(50.f, 50.f) + externalFieldDirection * 80.f, sf::Color::Blue)
        };

        window.draw(fieldLine, 2, sf::Lines);

        window.display();

        // std::cout << particles[0].angle << " "
        //   << particles[0].angularVelocity << " "
        //   << particles[0].angularAcceleration << "\n";

        //For radial distribution function
        // frameCounter++;
        // if (frameCounter % 30 == 0) {
        //     auto currentGR = computeRadialDistribution(particles, width, height, grBins, rMax);

        //     for (int i = 0; i < grBins; ++i) {
        //         grAverage[i] += currentGR[i];
        //         float value = grAverage[i] / grSamples;
        //     }

        //     grSamples++;

        //     if (frameCounter % 300 == 0 && !gr.empty()) {
        //         for (int i = 0; i < 10; ++i) {
        //             std::cout << i << " " << gr[i] << "\n";
        //         }
        //         std::cout << "----\n";
        //     }
        // }
    }
} 