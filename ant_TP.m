clc;
clear;
close all;

% Number of cities (nodes)
n = 50; 
coordinates = 100 * rand(n, 2); % Randomly generate city locations in 100x100 grid

% Distance matrix calculation
dist = zeros(n, n);
for i = 1:n
    for j = 1:n
        if i ~= j
            dist(i, j) = sqrt(sum((coordinates(i, :) - coordinates(j, :)).^2));
        else
            dist(i, j) = inf; % No self-loop
        end
    end
end

% ACO Parameters
num_ants = 20;         % Number of ants
alpha = 1;             % Influence of pheromone
beta = 1;              % Influence of heuristic (distance)
evaporation_rate = 0.5;% Pheromone evaporation rate
Q = 100;               % Pheromone deposit factor
num_iterations = 100;  % Number of iterations

% Initialize pheromone matrix
tau = ones(n, n); 

best_length = inf; % Initialize best tour length
best_tour = [];

% ACO Algorithm
for iter = 1:num_iterations
    tours = zeros(num_ants, n);  % Paths taken by ants
    tour_lengths = zeros(num_ants, 1);
    
    % Each ant constructs a tour
    for ant = 1:num_ants
        visited = false(1, n);
        start_node = randi(n);
        current_node = start_node;
        visited(current_node) = true;
        tour = [current_node];
        
        for step = 2:n
            % Probability selection for next city
            probs = (tau(current_node, :) .^ alpha) .* ((1 ./ dist(current_node, :)) .^ beta);
            probs(visited) = 0; % Remove visited cities
            probs = probs / sum(probs); % Normalize
            
            next_node = find(rand < cumsum(probs), 1);
            tour = [tour, next_node];
            visited(next_node) = true;
            current_node = next_node;
        end
        tour = [tour, start_node]; % Return to start
        
        % Compute tour length
        tour_length = sum(diag(dist(tour(1:end-1), tour(2:end))));
        
        % Store results
        tours(ant, :) = tour(1:end-1);
        tour_lengths(ant) = tour_length;
        
        % Update best solution
        if tour_length < best_length
            best_length = tour_length;
            best_tour = tour;
        end
    end
    
    % Pheromone update
    tau = (1 - evaporation_rate) * tau; % Evaporation
    
    for ant = 1:num_ants
        for j = 1:n
            a = tours(ant, j);
            b = tours(ant, mod(j, n) + 1);
            tau(a, b) = tau(a, b) + Q / tour_lengths(ant);
            tau(b, a) = tau(a, b); % Symmetric TSP
        end
    end
    
    % Display progress
    fprintf('Iteration %d, Best Tour Length: %.2f\n', iter, best_length);
end

% Plot results
figure;
hold on;
plot(coordinates(:, 1), coordinates(:, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
for i = 1:n
    text(coordinates(i, 1) + 2, coordinates(i, 2), sprintf('%d', i), 'FontSize', 10);
end
plot([coordinates(best_tour, 1); coordinates(best_tour(1), 1)], ...
     [coordinates(best_tour, 2); coordinates(best_tour(1), 2)], ...
     'b-', 'LineWidth', 2);
title(sprintf('Best Tour Length: %.2f', best_length));
grid on;
hold off;
