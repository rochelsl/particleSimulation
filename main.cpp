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

//For grid list, performance boost
float cellSize = rc;
int nx = static_cast<int>(width / cellSize);
int ny = static_cast<int>(height / cellSize);
std::vector<std::vector<std::vector<int>>> grid(nx, std::vector<std::vector<int>>(ny));

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;
    sf::Vector2f magneticMoment;
    float radius;
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

                            float r2 = r.x * r.x + r.y * r.y;
                            if (r2 > rc2 || r2 < 1e-6f) continue;

                            float sr2 = (sigma * sigma) / r2;
                            float sr6 = sr2 * sr2 * sr2;
                            float sr12 = sr6 * sr6;

                            // Scalar multiplying r-vector.
                            float f = 24.f * epsilon * (2.f * sr12 - sr6) / r2;
                            f -= f_shift;

                            sf::Vector2f force = r * f;

                            // mass = 1
                            p1.acceleration -= force;
                            p2.acceleration += force;
                        }
                    }
                }
            }
        }
    }
}

void integrate(std::vector<Particle>& particles, float dt) {
    for (auto& p : particles) {
        p.position += p.velocity * dt + 0.5f * p.acceleration * dt * dt;
        applyPBC(p.position, width, height);
        p.velocity += 0.5f * p.acceleration * dt;
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
    std::uniform_real_distribution<float> vel(-100.f, 100.f);

    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = 5.f;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        p.velocity = {vel(rng), vel(rng)};
        p.acceleration = {0.f, 0.f};
        p.magneticMoment = randomUnitVector(rng); //magnetization vector
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
        for (auto& p : particles)
            p.velocity += 0.5f * p.acceleration * dt;
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

        window.display();

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