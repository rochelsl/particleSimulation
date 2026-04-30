#include <SFML/Graphics.hpp>
#include "constants.hpp"
#include "mouse_field_controller.hpp"
#include "renderer.hpp"
#include "simulation.hpp"

int main() {
    sf::RenderWindow window(sf::VideoMode(sim::kWindowWidth, sim::kWindowHeight), "Dipolar Langevin Particle Simulation");
    window.setFramerateLimit(144);

    sim::Simulation simulation;
    sim::Renderer renderer;
    //sim::MouseFieldController fieldController;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            } else if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) {
                window.close();
            // } else if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
            //     fieldController.beginDrag(window.mapPixelToCoords({event.mouseButton.x, event.mouseButton.y}));
            // } else if (event.type == sf::Event::MouseMoved && fieldController.isDragging()) {
            //     fieldController.updateDrag(window.mapPixelToCoords({event.mouseMove.x, event.mouseMove.y}));
            // } else if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
            //     fieldController.endDrag(window.mapPixelToCoords({event.mouseButton.x, event.mouseButton.y}));
            }
        }
        //simulation.step(sim::kTimeStep, fieldController.fieldDirection(), fieldController.fieldStrength());
        //renderer.draw(window, simulation.particles(), fieldController);
    }
    return 0;
}