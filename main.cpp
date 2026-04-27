#include <SFML/Graphics.hpp> 
#include <vector>
#include <random>

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
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

int main() { 
    auto window = sf::RenderWindow{{1920u, 1080u}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    //Initialize many particles at random positions
    std::vector<Particle> particles;

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> posX(0.f, 1920.f);
    std::uniform_real_distribution<float> posY(0.f, 1080.f);
    std::uniform_real_distribution<float> vel(-100.f, 100.f);

    for (int i = 0; i < 2000; ++i) {
        Particle p;
        p.radius = 1.f;
        p.position = {posX(rng), posY(rng)};
        p.velocity = {vel(rng), vel(rng)};
        particles.push_back(p);
    }

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    sf::Clock clock;

    while (window.isOpen()) { 
        float dt = clock.restart().asSeconds();
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
        for (auto& p : particles) {
            p.position += p.velocity * dt;
            handleWallCollision(p, window.getSize());
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