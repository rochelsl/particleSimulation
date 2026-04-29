# (magnetic) particleSimulation

A simple simulation of particles interacting through a Lennard-Jones potential written in C++ and visualized using SFML (https://github.com/SFML/cmake-sfml-project).
Recent additions:
* magnetic dipole moments for each particles that is visualized by a red arrow
* dynamic, vector-based magnetic dipole-dipole interaction
* external magnetic background field
* the mouse can be used to control the direction and strength of the applied field

The simulation runs in real time and with 5000 particles at 30 fps on a laptop.

## Overview

Each particle is modeled as:

* A finite-size sphere (excluded volume)
* A periodic boundary conditions is applied
* A magnetic dipole with orientation
* Subject to thermal fluctuations (temperature)
* Evolving under translational and rotational dynamics

The simulation is visualized in real-time using SFML and allows interactive control of the external magnetic field.

## Controls
Left-click + drag → set external magnetic field (direction & strength)
Right-click → select particle (field visualization)
Toggle keys → enable/disable overlays (depending on version)

## Physical Model
### 1. Lennard-Jones Potential (Short-Range Interaction)

Used to model excluded volume and weak cohesion:

* Strong repulsion at short distances
* Weak attraction at intermediate distances

This prevents particles from collapsing into each other while allowing clustering.

### 2. Dipole-Dipole Interaction (Long-Range Anisotropic Force)

Each particle carries a magnetic moment. The interaction depends on relative orientation:

* Head-to-tail alignment → attraction → chain formation
* Side-by-side alignment → repulsion

### 3. External Magnetic Field

An external field:

* Applies a torque aligning dipoles
* Biases the system toward global alignment
* Promotes chain formation along the field direction

In the simulation:

Direction and strength can be set interactively (mouse drag)


### 4. Langevin Dynamics (Thermal Motion)

The system includes thermal fluctuations via a Langevin thermostat:

* Damping: slows motion (viscous medium)
* Noise: random kicks proportional to temperature

This models a particle suspended in a fluid (Brownian motion).

## Emergent Behavior

Depending on parameters, the system shows:

* Disordered gas-like motion (high temperature)
* Compact clusters (dominant LJ attraction)
* Dipolar chains (dominant dipole interaction + field)

A simulation of 5000 magnetic particles (LJ epsilon = 0.05, sigma = 8.0, particle radius = 0.5 * sigma, external field = 10000 in horizontal direction, simulation time step = 0.01, Langevin parameters: mass of each particle = 0.5, moment of inertia = 0.5, gamma translation = 2.0, gamma rotation = 8.0, initial temperature = 10.0, final temperature = 2.0, cooling time = 20.0)
<img width="2492" height="1406" alt="magnetic_particle_simulation_30_min_hor" src="https://github.com/user-attachments/assets/2f0ba8b9-6efc-4e91-8093-a3cd23a98e21" />

## Key Parameters
### Temperature (T)
* Controls thermal noise strength
* High → disorder
* Low → stable structures
### Dipole Strength (dipoleStrength)
* Strength of magnetic interaction
* High → strong chaining and alignment
* Low → LJ-dominated clustering
### External Field (H)
* Aligns dipoles globally
* Strong field → straight chains
* Weak field → random orientation
### Lennard-Jones Depth (epsilon)
* Controls isotropic attraction
* High → blob-like clusters
* Low → cleaner dipolar chains
### Damping (gammaT, gammaR)
* Controls how quickly motion relaxes
* High → overdamped (smooth, slow)
* Low → more inertial motion
### Time Step (dt)
* Numerical stability parameter
* Too large → overlaps / instability
* Smaller → more accurate, slower

## Dependencies
* C++17
* SFML
