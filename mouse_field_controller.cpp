#include "mouse_field_controller.hpp"

#include <algorithm>

#include "constants.hpp"
#include "math_utils.hpp"

namespace sim {

void MouseFieldController::beginDrag(const sf::Vector2f& position) {
    isDragging_ = true;
    dragStart_ = position;
    dragCurrent_ = position;
    updateFieldFromDrag();
}

void MouseFieldController::updateDrag(const sf::Vector2f& position) {
    dragCurrent_ = position;
    updateFieldFromDrag();
}

void MouseFieldController::endDrag(const sf::Vector2f& position) {
    dragCurrent_ = position;
    updateFieldFromDrag();
    isDragging_ = false;
}

bool MouseFieldController::isDragging() const { return isDragging_; }

sf::Vector2f MouseFieldController::dragStart() const { return dragStart_; }

sf::Vector2f MouseFieldController::dragCurrent() const { return dragCurrent_; }

sf::Vector2f MouseFieldController::fieldDirection() const { return normalize(fieldDirection_); }

float MouseFieldController::fieldStrength() const { return fieldStrength_; }

void MouseFieldController::updateFieldFromDrag() {
    const sf::Vector2f drag = dragCurrent_ - dragStart_;
    const float dragLength = length(drag);

    if (dragLength < kMinDragLengthForField) {
        fieldStrength_ = 0.f;
        return;
    }

    fieldDirection_ = drag / dragLength;
    fieldStrength_ = std::min(kMaxExternalFieldStrength, kFieldStrengthPerPixel * dragLength);
}

}  // namespace sim
