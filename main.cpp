#include <SFML/Graphics.hpp> 
#include <vector>
#include <random>
#include <cmath>

//Window size
const int width = 1920;
const int height = 1080;

//Initialize parameters for LJ-potential, cutoff and grid-calculations
const float epsilon = 1.0f;
const float sigma = 10.0f;
const float rc = 2.5f * sigma;
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
    float radius;
};

void handleWallCollision(Particle& p, const sf::Vector2u& size) {
    if (p.position.x - p.radius < 0.f) {
        p.position.x = p.radius;
        p.velocity.x *= -1.f;
    } else if (p.position.x + p.radius > size.x) {
        p.position.x = size.x - p.radius;
        p.velocity.x *= -1.f;
    }

    if (p.position.y - p.radius < 0.f) {
        p.position.y = p.radius;
        p.velocity.y *= -1.f;
    } else if (p.position.y + p.radius > size.y) {
        p.position.y = size.y - p.radius;
        p.velocity.y *= -1.f;
    }
}

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

int main() { 
    auto window = sf::RenderWindow{{width, height}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    std::vector<Particle> particles;

    const int N = 10000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-10.f, 10.f);

    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = 2.f;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        p.velocity = {vel(rng), vel(rng)};
        p.acceleration = {0.f, 0.f};
        particles.push_back(p);
    }

    // Initial acceleration for velocity Verlet.
    computeForces(particles);

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    //sf::Clock clock;
    //float dt = clock.restart().asSeconds();
    const float dt = 0.001f;

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
        }
        window.display();
    }
} 