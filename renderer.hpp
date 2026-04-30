#pragma once

#include <SFML/Graphics.hpp>

#include "mouse_field_controller.hpp"
#include "particle.hpp"

namespace sim {

class Renderer {
public:
    void draw(sf::RenderWindow& window, const std::vector<Particle>& particles, const MouseFieldController& fieldController);

private:
    sf::CircleShape shape_;
};

}  // namespace sim
