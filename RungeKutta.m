clear;
clc;

% Constants
G = 6.67430e-11;  % Gravitational constant [m^3 kg^-1 s^-2]
c = 299792458;    % Speed of light [m/s]
M = 4.0e6 * 1.989e30;  % Mass of the black hole (4 million solar masses)
a = 0.9 * c;      % Spin parameter of the black hole (close to maximal rotation)

% Initial conditions for a test particle (in spherical coordinates)
r0 = 3.0e10;      % Initial radial distance (in meters)
theta0 = pi/2;    % Initial angle in the equatorial plane (in radians)
phi0 = 0;         % Initial azimuthal angle (in radians)

v_r0 = 0;         % Initial radial velocity [m/s]
v_theta0 = 0;     % Initial angular velocity in theta direction [m/s]
v_phi0 = 2e4;     % Initial angular velocity in phi direction [m/s]

% Time parameters
t_max = 1e5;      % Maximum time for simulation [s]
dt = 1;           % Time step for integration [s]

% Initial position and velocity vectors
r = r0;
theta = theta0;
phi = phi0;

v_r = v_r0;
v_theta = v_theta0;
v_phi = v_phi0;

% Parameters derived from Kerr geodesics
E = (1 - 2*G*M/r + a^2/r^2) * v_r - (2*G*M/r) * a * sin(theta0)^2 * v_phi;  % Energy
L = r^2 * v_phi - a * (1 - 2*G*M/r + a^2/r^2) * v_r;  % Angular momentum

% Time vector
time = 0:dt:t_max;
orbits = zeros(length(time), 3);  % To store r, theta, and phi over time

% Runge-Kutta integration of the Kerr geodesics
for i = 1:length(time)
    % Store current values for plotting
    orbits(i, :) = [r, theta, phi];
    
    % Computing the effective potential for radial motion
    V_eff = (E^2 - (1 - 2*G*M/r + a^2/r^2) * (1 + L^2/r^2)) / (1 - 2*G*M/r + a^2/r^2);
    
    % Radial motion equation (effective potential)
    drdt = sqrt(abs(V_eff));  % Radial velocity
    
    % Runge-Kutta method for the integration of geodesic equations
    k1_r = dt * drdt;
    k1_phi = dt * (L/r^2);  % Angular velocity in phi direction

    % Updating r and phi using Runge-Kutta
    r_new = r + k1_r;
    phi_new = phi + k1_phi;

    % Updating the state variables
    r = r_new;
    phi = phi_new;
end

% Ploting the orbit showing the precession
figure;
plot3(orbits(:,1).*sin(orbits(:,2)), orbits(:,1).*cos(orbits(:,2)), orbits(:,3));
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Precessing Orbit Around a Rotating Black Hole');
grid on;


%% Stability Testing: Varying initial conditions for stability analysis

% Define multiple initial conditions for testing stability
initial_conditions = [3.0e10, 3.5e10, 4.0e10]; % Different initial radii
v_r_conditions = [0, 1e3, 2e3]; % Different radial velocities
v_phi_conditions = [2e4, 2.5e4, 3e4]; % Different angular velocities

% Time vector
time = 0:dt:t_max;

% Ploting for testing orbit stability
figure;
hold on;
for r_init = initial_conditions
    for v_r_init = v_r_conditions
        for v_phi_init = v_phi_conditions
            % Set initial conditions for each test
            r = r_init;
            phi = phi0;
            v_r = v_r_init;
            v_phi = v_phi_init;

            % Computing the initial energy and angular momentum
            E = (1 - 2*G*M/r + a^2/r^2) * v_r - (2*G*M/r) * a * sin(theta0)^2 * v_phi;
            L = r^2 * v_phi - a * (1 - 2*G*M/r + a^2/r^2) * v_r;

            % Iterating through the simulation using Runge-Kutta
            orbits = zeros(length(time), 3);
            for i = 1:length(time)
                % Storing current values for plotting
                orbits(i, :) = [r, theta, phi];
                
                % Radial motion equation
                V_eff = (E^2 - (1 - 2*G*M/r + a^2/r^2) * (1 + L^2/r^2)) / (1 - 2*G*M/r + a^2/r^2);
                drdt = sqrt(abs(V_eff));  % Radial velocity
                r = r + dt * drdt;  % Radial update
                
                % Updating phi
                phi = phi + dt * (L/r^2);  % Angular velocity in phi direction
            end
            
            % Ploting the orbit for the current set of initial conditions
            plot3(orbits(:,1).*sin(orbits(:,2)), orbits(:,1).*cos(orbits(:,2)), orbits(:,3));
        end
    end
end
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Orbit Stability Tests in Kerr Metric');
grid on;

