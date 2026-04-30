#pragma once

#include <SFML/Graphics.hpp>

namespace sim {

constexpr float PI = 3.14159265358979323846f;

constexpr int kWindowWidth = 1920;
constexpr int kWindowHeight = 1080;

constexpr float kEpsilon = 0.05f;
constexpr float kSigma = 8.0f;

constexpr float kParticleRadius = 0.5f * kSigma;
constexpr float kContactDistance = 3.0f * kParticleRadius;
constexpr float kDipoleSoftCoreDistance = kContactDistance;

constexpr float kMaxPairForce = 50000.0f;

constexpr float kLjCutoff = 2.5f * kSigma;
constexpr float kDipoleCutoff = 8.0f * kSigma;

constexpr float kDipoleStrength = 5000.0f;

constexpr float kFieldStrengthPerPixel = 20.0f;
constexpr float kMaxExternalFieldStrength = 10000.0f;
constexpr float kMinDragLengthForField = 2.0f;

constexpr float kMass = 0.5f;
constexpr float kMomentOfInertia = 0.5f;
constexpr float kGammaT = 2.0f;
constexpr float kGammaR = 8.0f;
constexpr float kInitialTemperature = 10.0f;
constexpr float kFinalTemperature = 2.0f;
constexpr float kCoolingTime = 20.0f;

constexpr int kParticleCount = 5000;
constexpr float kTimeStep = 0.01f;

}  // namespace sim
