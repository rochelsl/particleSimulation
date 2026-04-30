#pragma once

#include <SFML/Graphics.hpp>

namespace sim {

struct Particle {
    sf::Vector2f position{};
    sf::Vector2f velocity{};
    sf::Vector2f acceleration{};

    sf::Vector2f magneticMoment{};
    float angle = 0.f;
    float angularVelocity = 0.f;
    float angularAcceleration = 0.f;

    float radius = 5.f;
};

}  // namespace sim
