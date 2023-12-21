close all
clear 
clc

% Specify the file path
file_path = 'J:\MATLAB\NIN2272\New folder\TEST_1_L.txt';

% Read point cloud data from the file
data = readmatrix(file_path); % Assumes space-delimited values in the file

% Extract x, y, and z coordinates
x = data(:, 1)* 25.4;
y = data(:, 2)* 25.4;
z = data(:, 3)* 25.4;

% Calculate first and second derivatives
dt = 0.001; % Time step
vx = diff(x) / dt * 60; % First derivative of x (speed in x-direction)
vy = diff(y) / dt * 60; % First derivative of y (speed in y-direction)
vz = diff(z) / dt * 60; % First derivative of z (speed in z-direction)

% unit conversion before taking the 2nd derivative
% vx = vx / 1000 * 60;
% vy = vy / 1000 * 60;
% vx = vz / 1000 * 60;

% Exclude the last point for the second derivative calculation
ax = diff(vx) / dt; % Second derivative of x (acceleration in x-direction)
ay = diff(vy) / dt; % Second derivative of y (acceleration in y-direction)
az = diff(vz) / dt; % Second derivative of z (acceleration in z-direction)

% Time vector (excluding the last two points due to differentiation)
t = (0:dt:(length(x)-2)*dt);

% Plot the results
figure;
subplot(2, 1, 1);
plot(t, sqrt(vx.^2 + vy.^2 + vz.^2), 'b-', 'LineWidth', 1.5);
title('Feedrate vs Time');
xlabel('Time (s)');
ylabel('Speed (mm/min)');

% Time vector for acceleration (excluding the last point)
t_acc = (0:dt:(length(x)-3)*dt);

subplot(2, 1, 2);
plot(t_acc, sqrt(ax.^2 + ay.^2 + az.^2), 'r-', 'LineWidth', 1.5);
title('Acceleration vs Time');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

