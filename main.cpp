#include <SFML/Graphics.hpp> 
#include <vector>
#include <random>

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

void computeForces(std::vector<Particle>& particles) {
    const float epsilon = 1.0f;
    const float sigma = 10.0f;
    const float rc = 2.5f * sigma;
    const float rc2 = rc * rc;

    float inv_rc2 = 1.f / rc2;
    float inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
    float inv_rc12 = inv_rc6 * inv_rc6;

    // force at cutoff
    float f_shift = 24.f * epsilon * inv_rc2 * (2.f * inv_rc12 - inv_rc6);

    // reset accelerations
    for (auto& p : particles)
        p.acceleration = {0.f, 0.f};

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            auto& p1 = particles[i];
            auto& p2 = particles[j];

            sf::Vector2f r = p2.position - p1.position;
            float r2 = r.x*r.x + r.y*r.y;

            if (r2 > rc2 || r2 < 1e-6f)
                continue;

            float inv_r2 = 1.f / r2;
            float inv_r6 = inv_r2 * inv_r2 * inv_r2;
            float inv_r12 = inv_r6 * inv_r6;

            float f = 24.f * epsilon * inv_r2 * (2.f * inv_r12 - inv_r6);

            // subtract cutoff force
            f -= f_shift;

            sf::Vector2f force = r * f;

            // Newton's 3rd law
            p1.acceleration -= force;
            p2.acceleration += force;
        }
    }
}

void integrate(std::vector<Particle>& particles, float dt) {
    for (auto& p : particles) {
        p.position += p.velocity * dt + 0.5f * p.acceleration * dt * dt;
        p.velocity += 0.5f * p.acceleration * dt;
    }
}

int main() { 
    auto window = sf::RenderWindow{{1920u, 1080u}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    //Initialize many particles at random positions
    std::vector<Particle> particles;

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> posX(50.f, 1880.f);
    std::uniform_real_distribution<float> posY(50.f, 1030.f);
    std::uniform_real_distribution<float> vel(-100.f, 100.f);

    for (int i = 0; i < 1000; ++i) {
        Particle p;
        p.radius = 3.f;
        p.position = {posX(rng), posY(rng)};
        p.velocity = {vel(rng), vel(rng)};
        particles.push_back(p);
    }

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    sf::Clock clock;

    while (window.isOpen()) { 
        //float dt = clock.restart().asSeconds();
        const float dt = 0.005f;
        // 1. integrate positions + half velocity
        integrate(particles, dt);
        // 2. recompute forces
        computeForces(particles);
        // 3. finish velocity update
        for (auto& p : particles)
            p.velocity += 0.5f * p.acceleration * dt;
        // 4. walls (still needed)
        for (auto& p : particles)
            handleWallCollision(p, window.getSize());

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

        window.clear();
        //Collision with the wall check
        for (const auto& p : particles) {
            shape.setRadius(p.radius);
            shape.setOrigin({p.radius, p.radius});
            shape.setPosition(p.position);
            window.draw(shape);
        }

        //Collision with each other check
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                auto& p1 = particles[i];
                auto& p2 = particles[j];

                sf::Vector2f delta = p2.position - p1.position;
                float dist2 = delta.x * delta.x + delta.y * delta.y;
                float minDist = p1.radius + p2.radius;
                if (dist2 < 1e-8f) continue;
                if (dist2 < minDist * minDist) {
                    float dist = std::sqrt(dist2);
                    sf::Vector2f normal = delta / dist;

                    // Overlap correction
                    float percent = 0.8f;   // correction strength
                    float slop = 0.01f;     // tolerance

                    float overlap = std::max(minDist - dist - slop, 0.f);
                    sf::Vector2f correction = normal * (percent * overlap * 0.5f);
                    p1.position -= correction;
                    p2.position += correction;

                    // Relative velocity
                    sf::Vector2f relVel = p2.velocity - p1.velocity;
                    float velAlongNormal = relVel.x * normal.x + relVel.y * normal.y;

                    //moving apart -> skip
                    if (velAlongNormal > 0.f)
                        continue;

                    // Elastic response
                    float j = -2.f * velAlongNormal / 2.f; // simplifies to -velAlongNormal

                    sf::Vector2f impulse = j * normal;

                    p1.velocity -= impulse;
                    p2.velocity += impulse;


                }
            }
        }
        window.display();
    }
} 