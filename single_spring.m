close all;
clear;

%% Define global parameters
g = 10;  % Acceleration due to gravity (m/s^2)
m = 1.0;  % Mass of the particle (kg)
L0 = 1.0;  % Rest length of the spring (meters)
initial_theta = 0.5;  % Initial angle (radians)
initial_L = 0;  % Initial extension/compression (meters)

% Different combinations of spring constants (k) and damping constants (beta)
params = [3.0, 0.6; 6.0, 1.0; 10.0, 2.0; 8.0, 4.0];

%% Simulation setup
simOut = cell(size(params, 1), 1);
for i = 1:size(params, 1)
    k = params(i, 1);  % Spring constant
    beta = params(i, 2);  % Damping constant

    % Assign parameters to workspace for Simulink
    assignin('base', 'k', k);
    assignin('base', 'beta', beta);
    assignin('base', 'g', g);
    assignin('base', 'm', m);
    assignin('base', 'L0', L0);
    assignin('base', 'initial_theta', initial_theta);
    assignin('base', 'initial_L', initial_L);

    % Run Simulink model
    simOut{i} = sim('spring_simulator', 'SimulationMode', 'normal');
end

%% Prepare figure for animation
h_fig = figure;
hold on;
grid on;
xlabel('X Position (meters)');
ylabel('Y Position (meters)');
title('Spring Mass Systems with Different Parameters');
% axis equal;

% Plot references and prepare dynamic objects
ref_x = 0.0;
ref_y = 0.0;
plot(ref_x, ref_y, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');  % Reference point

% Custom colors for each simulation
colors = [1 0 0; 0 0 1; 1 0.5 0; 0 1 0];  % Red, Blue, Orange, Green
h_plot = gobjects(size(params, 1), 1);
h_line = gobjects(size(params, 1), 1);
legends = strings(size(params, 1), 1);

x_all = [];
y_all = [];

for i = 1:size(params, 1)
    time = simOut{i}.position.time;
    position_values = simOut{i}.position.signals.values;
    x = position_values(:, 1);
    y = position_values(:, 2);

    x_all = [x_all; x];  % Collect all x values
    y_all = [y_all; y];  % Collect all y values

    h_plot(i) = plot(x(1), y(1), 'o', 'MarkerSize', 6, 'MarkerFaceColor', colors(i, :));
    h_line(i) = line([x(1), ref_x], [y(1), ref_y], 'Color', colors(i, :), 'LineStyle', '-');
    legends(i) = sprintf('k = %.1f, beta = %.1f', params(i,1), params(i,2));
end

legend(h_plot, legends, 'Location', 'best');

% Set fixed axis limits
xlim([-4, 4]);
ylim([-5, 0]);

%% Animate all systems in the same plot
maxTimeSteps = max(cellfun(@(c) numel(c.position.time), simOut));  % Find the longest time vector
frameSkip = 5;  % Update the plot every 'frameSkip' frames for faster animation
for k = 1:frameSkip:maxTimeSteps
    for i = 1:size(params, 1)
        if k <= length(simOut{i}.position.time)
            position_values = simOut{i}.position.signals.values;
            x = position_values(:, 1);
            y = position_values(:, 2);

            % Update the position of markers and lines
            set(h_plot(i), 'XData', x(k), 'YData', y(k));
            set(h_line(i), 'XData', [x(k), ref_x], 'YData', [y(k), ref_y]);
        end
    end
    drawnow;
end
