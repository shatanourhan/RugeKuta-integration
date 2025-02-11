# Black Hole Orbit Simulation in Kerr Metric
Overview
This MATLAB script simulates the orbits of test particles (e.g., stars or planets) around a rotating black hole using the Kerr metric. The simulation incorporates Runge-Kutta numerical integration to solve the geodesic equations, allowing the study of orbital precession and stability.

Features
Uses the Kerr metric to model the gravitational field of a rotating black hole.

Implements numerical integration via the Runge-Kutta method.

Simulates orbital trajectories and visualizes the precession effects.

Tests the stability of different orbits by varying initial conditions.

3D visualization of orbits in the black hole's gravitational field.

Requirements
MATLAB (R2020a or later recommended)

Basic understanding of general relativity and geodesic equations

Installation
Download the black_hole_orbits.m MATLAB script.

Open MATLAB and navigate to the folder containing the script.

Run the script by typing:

black_hole_orbits
Parameters
The following key parameters can be adjusted:

M: Mass of the black hole (default is 4 million solar masses, similar to Sagittarius A*).

a: Spin parameter of the black hole (default is 0.9, indicating near-maximal rotation).

r0: Initial radial distance of the test particle.

v_r0, v_phi0: Initial radial and angular velocities.

t_max, dt: Maximum simulation time and time step for integration.

Output
A 3D plot showing the trajectory of the orbiting object.

A stability test plot showing variations in orbits based on different initial conditions.

Insights into orbital precession and effective potential variations in a Kerr spacetime.

Usage
Modify the initial conditions in the script to analyze different types of orbits:

r0 = 3.0e10;  % Change initial radial distance

v_phi0 = 2.5e4;  % Modify angular velocity
Run the script again to observe the impact of these changes on the orbit.

Future Improvements
Implementing adaptive step-size Runge-Kutta for improved accuracy.

Extending to include the effect of gravitational waves.

Simulating orbits in an extreme mass-ratio inspiral (EMRI) scenario.

Author
Developed by Shata Nourhan. For questions or contributions, feel free to contact me.

License
This project is open-source and distributed under the MIT License.
