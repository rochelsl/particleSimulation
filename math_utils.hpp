#pragma once

#include <SFML/Graphics.hpp>

namespace sim {

float dot(sf::Vector2f a, sf::Vector2f b);
float length(sf::Vector2f v);
float cross2D(sf::Vector2f a, sf::Vector2f b);
sf::Vector2f normalize(sf::Vector2f v);
sf::Vector2f neg(sf::Vector2f v);
sf::Vector2f momentFromAngle(float angle);

void applyPBC(sf::Vector2f& pos, float width, float height);
sf::Vector2f minimumImage(sf::Vector2f r, float width, float height);

}  // namespace sim
