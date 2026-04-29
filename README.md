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

# Mathematical model
## Equations of motion
### Translational motion

%m \frac{d\mathbf{v}_i}{dt} = \mathbf{F}_i - \gamma_T \mathbf{v}_i + \sqrt{2 \gamma_T k_B T m}, \boldsymbol{\eta}_i(t)%

### Rotational motion

%I \frac{d\omega_i}{dt} = \tau_i - \gamma_R \omega_i + \sqrt{2 \gamma_R k_B T I}, \eta_i(t)%

* ( $m$ ): particle mass
* ( $I$ ): moment of inertia
* ( $\gamma_T$, $\gamma_R$ ): translational and rotational damping
* ( $T$ ): temperature
* ( $\boldsymbol{\eta}(t)$ ): Gaussian white noise

## 2. Lennard-Jones Interaction

Pairwise isotropic interaction:

$U_{\mathrm{LJ}}(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]$

Force:
$[
\mathbf{F}*{ij}^{\mathrm{LJ}} = -\nabla U*{\mathrm{LJ}}(r_{ij})
]$

Parameters:

* ( $\epsilon$ ): interaction strength
* ( $\sigma$ ): particle diameter scale

## 3. Dipole–Dipole Interaction

Each particle carries a dipole moment ( \mathbf{m}_i ).

### Interaction energy

$[
U_{ij}^{\mathrm{dip}} = \frac{\mu_0}{4\pi r^3}
\left[
\mathbf{m}_i \cdot \mathbf{m}_j

* 3 (\mathbf{m}_i \cdot \hat{\mathbf{r}})
  (\mathbf{m}_j \cdot \hat{\mathbf{r}})
  \right]
  ]$

### Force

$[
\mathbf{F}*{ij}^{\mathrm{dip}} = -\nabla U*{ij}^{\mathrm{dip}}
]$

### Torque

$[
\boldsymbol{\tau}_i^{\mathrm{dip}} = \mathbf{m}_i \times \mathbf{B}_j
]$

where ( $\mathbf{B}_j$ ) is the magnetic field generated by dipole ( $j$ ).

## 4. Magnetic Field of a Dipole

$\mathbf{B}(\mathbf{r}) = \frac{\mu_0}{4\pi r^3} \left[ 3(\mathbf{m} \cdot \hat{\mathbf{r}}) \hat{\mathbf{r}} - \mathbf{m} \right]$

This field defines both:

* dipole–dipole forces
* torque alignment

## 5. External Magnetic Field

Uniform external field ( \mathbf{H} ):

### Torque on particle

$[
\boldsymbol{\tau}_i^{\mathrm{ext}} = \mathbf{m}_i \times \mathbf{H}
]$

### Energy

$[
U_i^{\mathrm{ext}} = -\mathbf{m}_i \cdot \mathbf{H}
]$

Effect:

* Aligns dipoles globally
* Drives chain formation

## 6. Total Force and Torque

$[
\mathbf{F}*i =
\sum*{j \ne i}
\left(
\mathbf{F}*{ij}^{\mathrm{LJ}} +
\mathbf{F}*{ij}^{\mathrm{dip}}
\right)
]

[
\boldsymbol{\tau}_i =
\boldsymbol{\tau}_i^{\mathrm{dip}} +
\boldsymbol{\tau}_i^{\mathrm{ext}}
]$

## 7. Dimensionless Control Ratios

Key behavior is governed by:

* Thermal vs dipolar energy:
 $ [
  \frac{U_{\mathrm{dip}}}{k_B T}
  ]$

* Dipolar vs isotropic interaction:
 $ [
  \frac{\text{dipoleStrength}}{\epsilon}
  ]$

These ratios determine whether the system forms:

* gas-like states
* isotropic clusters
* dipolar chains
