#include "renderer.hpp"

#include "constants.hpp"

namespace sim {

void Renderer::draw(sf::RenderWindow& window, const std::vector<Particle>& particles, const MouseFieldController& fieldController) {
    shape_.setFillColor(sf::Color::White);

    window.clear();

    for (const auto& p : particles) {
        shape_.setRadius(p.radius);
        shape_.setOrigin(p.radius, p.radius);
        shape_.setPosition(p.position);
        window.draw(shape_);

        const sf::Vertex line[] = {
            sf::Vertex(p.position, sf::Color::Red),
            sf::Vertex(p.position + p.magneticMoment * p.radius * 2.0f, sf::Color::Red)};
        window.draw(line, 2, sf::Lines);
    }

    const sf::Vector2f origin(50.f, 50.f);
    const float indicatorLength = 120.f * fieldController.fieldStrength() / kMaxExternalFieldStrength;
    const sf::Vertex fieldLine[] = {
        sf::Vertex(origin, sf::Color::Yellow),
        sf::Vertex(origin + fieldController.fieldDirection() * indicatorLength, sf::Color::Yellow)};
    window.draw(fieldLine, 2, sf::Lines);

    if (fieldController.isDragging()) {
        const sf::Vertex dragLine[] = {
            sf::Vertex(fieldController.dragStart(), sf::Color::Green),
            sf::Vertex(fieldController.dragCurrent(), sf::Color::Green)};
        window.draw(dragLine, 2, sf::Lines);
    }

    window.display();
}

}  // namespace sim
