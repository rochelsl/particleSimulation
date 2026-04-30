#pragma once

#include <SFML/Graphics.hpp>

namespace sim {

class MouseFieldController {
public:
    void beginDrag(const sf::Vector2f& position);
    void updateDrag(const sf::Vector2f& position);
    void endDrag(const sf::Vector2f& position);

    bool isDragging() const;
    sf::Vector2f dragStart() const;
    sf::Vector2f dragCurrent() const;

    sf::Vector2f fieldDirection() const;
    float fieldStrength() const;

private:
    void updateFieldFromDrag();

    bool isDragging_ = false;
    sf::Vector2f dragStart_{0.f, 0.f};
    sf::Vector2f dragCurrent_{0.f, 0.f};

    sf::Vector2f fieldDirection_{1.f, 0.2f};
    float fieldStrength_ = 100.0f;
};

}  // namespace sim
