#include "math_utils.hpp"

#include <cmath>

namespace sim {

float dot(sf::Vector2f a, sf::Vector2f b) { return a.x * b.x + a.y * b.y; }

float length(sf::Vector2f v) { return std::sqrt(dot(v, v)); }

float cross2D(sf::Vector2f a, sf::Vector2f b) { return a.x * b.y - a.y * b.x; }

sf::Vector2f normalize(sf::Vector2f v) {
    const float len = std::sqrt(dot(v, v));
    if (len < 1e-8f) return {0.f, 0.f};
    return v / len;
}

sf::Vector2f neg(sf::Vector2f v) { return {-v.x, -v.y}; }

sf::Vector2f momentFromAngle(float angle) { return {std::cos(angle), std::sin(angle)}; }

void applyPBC(sf::Vector2f& pos, float width, float height) {
    pos.x = std::fmod(pos.x, width);
    pos.y = std::fmod(pos.y, height);
    if (pos.x < 0.f) pos.x += width;
    if (pos.y < 0.f) pos.y += height;
}

sf::Vector2f minimumImage(sf::Vector2f r, float width, float height) {
    if (r.x > 0.5f * width) r.x -= width;
    if (r.x < -0.5f * width) r.x += width;
    if (r.y > 0.5f * height) r.y -= height;
    if (r.y < -0.5f * height) r.y += height;
    return r;
}

}  // namespace sim
