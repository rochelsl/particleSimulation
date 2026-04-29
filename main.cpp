#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

constexpr float PI = 3.14159265358979323846f;

// Window size
constexpr int width  = 1920;
constexpr int height = 1080;

// Lennard-Jones parameters
constexpr float epsilon = 1.0f;
constexpr float sigma   = 10.0f;
constexpr float ljCutoff = 2.5f * sigma;      // 1.0*sigma removes nearly all LJ attraction
constexpr float ljCutoff2 = ljCutoff * ljCutoff;

// Dipolar interaction parameters
constexpr float dipoleStrength = 100.0f;
constexpr float dipoleCutoff = 6.0f * sigma;  // finite cutoff for the otherwise long-ranged dipole force
constexpr float dipoleCutoff2 = dipoleCutoff * dipoleCutoff;
constexpr float rotationalDamping = 0.999f;

sf::Vector2f normalize(sf::Vector2f v) {
    const float len = std::sqrt(v.x * v.x + v.y * v.y);
    if (len < 1e-8f) return {0.f, 0.f};
    return v / len;
}

// Background magnetic field
const sf::Vector2f externalFieldDirection = normalize({1.0f, 0.2f});
constexpr float externalFieldStrength = 100.0f;

// Cell-list parameters. The 3x3 neighbor search is valid when cellSize >= largest cutoff.
const float cellSize = std::max(ljCutoff, dipoleCutoff);
const int nx = std::max(1, static_cast<int>(std::ceil(width  / cellSize)));
const int ny = std::max(1, static_cast<int>(std::ceil(height / cellSize)));
std::vector<std::vector<std::vector<int>>> grid(nx, std::vector<std::vector<int>>(ny));

struct Particle {
    sf::Vector2f position{};
    sf::Vector2f velocity{};
    sf::Vector2f acceleration{};

    // Magnetic moment orientation
    sf::Vector2f magneticMoment{};
    float angle = 0.f;
    float angularVelocity = 0.f;
    float angularAcceleration = 0.f;

    float radius = 5.f;
};

float dot(sf::Vector2f a, sf::Vector2f b) {
    return a.x * b.x + a.y * b.y;
}

float cross2D(sf::Vector2f a, sf::Vector2f b) {
    return a.x * b.y - a.y * b.x;
}

sf::Vector2f momentFromAngle(float angle) {
    return {std::cos(angle), std::sin(angle)};
}

sf::Vector2f randomUnitVector(std::mt19937& rng) {
    std::uniform_real_distribution<float> angleDist(0.f, 2.f * PI);
    return momentFromAngle(angleDist(rng));
}

void applyPBC(sf::Vector2f& pos, float w, float h) {
    pos.x = std::fmod(pos.x, w);
    pos.y = std::fmod(pos.y, h);
    if (pos.x < 0.f) pos.x += w;
    if (pos.y < 0.f) pos.y += h;
}

sf::Vector2f minimumImage(sf::Vector2f r, float w, float h) {
    if (r.x >  0.5f * w) r.x -= w;
    if (r.x < -0.5f * w) r.x += w;
    if (r.y >  0.5f * h) r.y -= h;
    if (r.y < -0.5f * h) r.y += h;
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

sf::Vector2f computeDipoleForceOnP2(
    sf::Vector2f r,       // vector from particle 1 to particle 2
    sf::Vector2f mu1,
    sf::Vector2f mu2,
    float strength
) {
    const float r2 = dot(r, r);
    if (r2 < 1e-6f) return {0.f, 0.f};

    const float dist = std::sqrt(r2);
    const sf::Vector2f rHat = r / dist;

    const float mu1_dot_mu2 = dot(mu1, mu2);
    const float mu1_dot_r   = dot(mu1, rHat);
    const float mu2_dot_r   = dot(mu2, rHat);

    const float prefactor = 3.f * strength / (dist * dist * dist * dist);

    const sf::Vector2f bracket =
        mu1_dot_r * mu2 +
        mu2_dot_r * mu1 +
        mu1_dot_mu2 * rHat -
        5.f * mu1_dot_r * mu2_dot_r * rHat;

    return prefactor * bracket;
}

sf::Vector2f dipoleField(sf::Vector2f r, sf::Vector2f mu, float strength) {
    const float r2 = dot(r, r);
    if (r2 < 1e-6f) return {0.f, 0.f};

    const float dist = std::sqrt(r2);
    const sf::Vector2f rHat = r / dist;
    const float muDotR = dot(mu, rHat);

    return strength * (3.f * muDotR * rHat - mu) / (dist * dist * dist);
}

void computeForces(std::vector<Particle>& particles) {
    for (auto& col : grid) {
        for (auto& cell : col) cell.clear();
    }

    const sf::Vector2f B_ext = externalFieldStrength * externalFieldDirection;

    for (auto& p : particles) {
        applyPBC(p.position, static_cast<float>(width), static_cast<float>(height));
        p.acceleration = {0.f, 0.f};

        // Include the external field in the same force/torque evaluation as the pair torques.
        p.angularAcceleration = cross2D(p.magneticMoment, B_ext);
    }

    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        const int cx = cellIndex(particles[i].position.x, nx, static_cast<float>(width));
        const int cy = cellIndex(particles[i].position.y, ny, static_cast<float>(height));
        grid[cx][cy].push_back(i);
    }

    // Force-shifted LJ scalar. The force vector is r * [scalar - scalar(rc)].
    const float sr2c = (sigma * sigma) / ljCutoff2;
    const float sr6c = sr2c * sr2c * sr2c;
    const float sr12c = sr6c * sr6c;
    const float ljScalarAtCutoff = 24.f * epsilon * (2.f * sr12c - sr6c) / ljCutoff2;

    for (int cx = 0; cx < nx; ++cx) {
        for (int cy = 0; cy < ny; ++cy) {
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    const int cx2 = (cx + dx + nx) % nx;
                    const int cy2 = (cy + dy + ny) % ny;

                    for (int i : grid[cx][cy]) {
                        for (int j : grid[cx2][cy2]) {
                            if (i >= j) continue;

                            Particle& p1 = particles[i];
                            Particle& p2 = particles[j];

                            sf::Vector2f r = p2.position - p1.position;
                            r = minimumImage(r, static_cast<float>(width), static_cast<float>(height));

                            const float r2 = dot(r, r);
                            if (r2 < 1e-6f) continue;

                            sf::Vector2f totalForce{0.f, 0.f};

                            if (r2 < ljCutoff2) {
                                const float inv_r2 = 1.f / r2;
                                const float sr2 = (sigma * sigma) * inv_r2;
                                const float sr6 = sr2 * sr2 * sr2;
                                const float sr12 = sr6 * sr6;

                                const float ljScalar =
                                    24.f * epsilon * inv_r2 * (2.f * sr12 - sr6)
                                    - ljScalarAtCutoff;

                                totalForce += r * ljScalar;
                            }

                            if (r2 < dipoleCutoff2) {
                                // Pair force
                                totalForce += computeDipoleForceOnP2(
                                    r,
                                    p1.magneticMoment,
                                    p2.magneticMoment,
                                    dipoleStrength
                                );

                                // Pair torques
                                const sf::Vector2f B_on_1 = dipoleField(-r, p2.magneticMoment, dipoleStrength);
                                const sf::Vector2f B_on_2 = dipoleField( r, p1.magneticMoment, dipoleStrength);

                                p1.angularAcceleration += cross2D(p1.magneticMoment, B_on_1);
                                p2.angularAcceleration += cross2D(p2.magneticMoment, B_on_2);
                            }

                            // Unit mass is assumed, so force == acceleration.
                            p1.acceleration -= totalForce;
                            p2.acceleration += totalForce;
                        }
                    }
                }
            }
        }
    }
}

void integrateFirstHalf(std::vector<Particle>& particles, float dt) {
    for (auto& p : particles) {
        p.position += p.velocity * dt + p.acceleration * (0.5f * dt * dt);
        p.velocity += p.acceleration * (0.5f * dt);

        p.angle += p.angularVelocity * dt + 0.5f * p.angularAcceleration * dt * dt;
        p.angularVelocity += 0.5f * p.angularAcceleration * dt;

        p.magneticMoment = momentFromAngle(p.angle);
        applyPBC(p.position, static_cast<float>(width), static_cast<float>(height));
    }
}

void integrateSecondHalf(std::vector<Particle>& particles, float dt) {
    for (auto& p : particles) {
        p.velocity += p.acceleration * (0.5f * dt);
        p.angularVelocity += 0.5f * p.angularAcceleration * dt;
        p.angularVelocity *= rotationalDamping;
        p.magneticMoment = momentFromAngle(p.angle);
    }
}

std::vector<float> computeRadialDistribution(
    const std::vector<Particle>& particles,
    float boxWidth,
    float boxHeight,
    int bins,
    float rMax
) {
    std::vector<float> hist(bins, 0.f);

    const int N = static_cast<int>(particles.size());
    const float dr = rMax / bins;
    const float area = boxWidth * boxHeight;
    const float rho = static_cast<float>(N) / area;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            sf::Vector2f r = particles[j].position - particles[i].position;
            r = minimumImage(r, boxWidth, boxHeight);

            const float dist = std::sqrt(dot(r, r));
            if (dist < rMax) {
                const int bin = static_cast<int>(dist / dr);
                if (bin >= 0 && bin < bins) hist[bin] += 2.f;
            }
        }
    }

    std::vector<float> g(bins, 0.f);
    for (int k = 0; k < bins; ++k) {
        const float rInner = k * dr;
        const float rOuter = rInner + dr;
        const float shellArea = PI * (rOuter * rOuter - rInner * rInner);
        const float idealCount = static_cast<float>(N) * rho * shellArea;
        if (idealCount > 0.f) g[k] = hist[k] / idealCount;
    }

    return g;
}

int main() {
    sf::RenderWindow window(sf::VideoMode(width, height), "Particle Simulation");
    window.setFramerateLimit(144);

    std::vector<Particle> particles;
    constexpr int N = 1000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-100.f, 100.f);
    std::uniform_real_distribution<float> angleDist(0.f, 2.f * PI);

    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    const int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    const int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    const float dx = static_cast<float>(width) / cols;
    const float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        const int ix = i % cols;
        const int iy = i / cols;

        Particle p;
        p.radius = 5.f;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        p.velocity = {vel(rng), vel(rng)};
        p.angle = angleDist(rng);
        p.angularVelocity = 10.f;
        p.magneticMoment = momentFromAngle(p.angle);

        particles.push_back(p);
    }

    // Initial acceleration for velocity Verlet.
    computeForces(particles);

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    constexpr float dt = 0.001f;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            } else if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) {
                window.close();
            }
        }

        integrateFirstHalf(particles, dt);
        computeForces(particles);
        integrateSecondHalf(particles, dt);

        window.clear();

        for (const auto& p : particles) {
            shape.setRadius(p.radius);
            shape.setOrigin(p.radius, p.radius);
            shape.setPosition(p.position);
            window.draw(shape);

            const sf::Vertex line[] = {
                sf::Vertex(p.position, sf::Color::Red),
                sf::Vertex(p.position + p.magneticMoment * p.radius * 2.0f, sf::Color::Red)
            };
            window.draw(line, 2, sf::Lines);
        }

        const sf::Vector2f origin(50.f, 50.f);
        const sf::Vertex fieldLine[] = {
            sf::Vertex(origin, sf::Color::Blue),
            sf::Vertex(origin + externalFieldDirection * 80.f, sf::Color::Blue)
        };
        window.draw(fieldLine, 2, sf::Lines);

        window.display();
    }

    return 0;
}