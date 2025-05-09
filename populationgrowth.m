% Data provided
year = [1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010];
population = [1750, 1860, 2070, 2300, 2520, 3020, 3700, 4440, 5270, 5980, 6790];

% Start time
tic;

% Classical exponential growth model
classical_model = @(P, t) population(1) * exp(P * (t - year(1)));

% Mittag-Leffler function (approximation)
mlf = @(alpha, z) arrayfun(@(z) sum(z.^(0:100) ./ gamma(alpha * (0:100) + 1)), z);

% Fractional model
fractional_model = @(params, t) population(1) * mlf(params(1), params(2) * (t - year(1)).^params(1));

% Different initial guesses
initial_guesses_classical = [0.005, 0.01, 0.015];
initial_guesses_fractional = [1.2, 0.003; 1.3, 0.002; 1.5 ,0.01];


t_fit = linspace(min(year), max(year), 1000);
figure;
hold on;
plot(year, population, '-x', 'DisplayName', 'Data');


% Preallocate storage
errors_classical = zeros(length(initial_guesses_classical), 1);
errors_fractional = zeros(size(initial_guesses_fractional, 1), 1);
optimized_P_classical = zeros(length(initial_guesses_classical), 1);
optimized_params_fractional = zeros(size(initial_guesses_fractional, 1), 2);


% Loop over different initial guesses for classical model
for i = 1:length(initial_guesses_classical)
    P0_classical = initial_guesses_classical(i);
    P_classical = lsqcurvefit(classical_model, P0_classical, year, population);
    N_classical = classical_model(P_classical, t_fit);
    
    % Store optimized parameter and error
    optimized_P_classical(i) = P_classical;
    errors_classical(i) = sum((population - classical_model(P_classical, year)).^2);

    % Plot classical model
    plot(t_fit, N_classical, 'DisplayName', sprintf('Classical (P0=%.3f, P=%.6f)', P0_classical, P_classical));
end


% Loop over different initial guesses for fractional model
for j = 1:size(initial_guesses_fractional, 1)
    params0_fractional = initial_guesses_fractional(j, :);
    params_fractional = lsqcurvefit(fractional_model, params0_fractional, year, population);
    N_fractional = fractional_model(params_fractional, t_fit);
    
    % Store optimized parameters and error
    optimized_params_fractional(j, :) = params_fractional;
    errors_fractional(j) = sum((population - fractional_model(params_fractional, year)).^2);

    % Plot fractional model
    plot(t_fit, N_fractional, '--', 'DisplayName', ...
        sprintf('Frac (α0=%.1f, P0=%.3f, α=%.3f, P=%.6f)', ...
        params0_fractional(1), params0_fractional(2), params_fractional(1), params_fractional(2)));
end

xlabel('Year');
ylabel('Population');
legend;
title('Comparison of Classical and Fractional Growth Models');
grid on;
hold off;


inset_axes = axes('Position', [0.64 0.2 0.25 0.25]); % [x, y, width, height]
box on; hold on; grid on;

% Define zoomed-in limits
xlim_range = [1927.04275 1927.04295];
ylim_range = [2016.39874 2016.39876];

% Plot data points inside the inset
plot(inset_axes, year, population, '-x', 'DisplayName', 'Data');

% Replot all models inside the zoomed-in inset
for i = 1:length(initial_guesses_classical)
    plot(inset_axes, t_fit, classical_model(optimized_P_classical(i), t_fit), 'LineWidth', 1);
end
for j = 1:size(initial_guesses_fractional, 1)
    plot(inset_axes, t_fit, fractional_model(optimized_params_fractional(j, :), t_fit), '--', 'LineWidth', 1);
end

% Adjust inset plot properties
xlim(inset_axes, xlim_range);
ylim(inset_axes, ylim_range);
title(inset_axes, 'Zoomed View');
xlabel(inset_axes, 'Year');
ylabel(inset_axes, 'Population');


% Stop time
stop_time = toc;
fprintf('Execution time: %f seconds\n', stop_time);

fprintf('\nErrors and Efficiency Gains:\n\n');


for i = 1:length(initial_guesses_classical)
    fprintf('Classical (P0=%.3f) Error: %f, Optimized P: %.9f\n', ...
        initial_guesses_classical(i), errors_classical(i), optimized_P_classical(i));
end

for j = 1:size(initial_guesses_fractional, 1)
    fprintf('Fractional (α0=%.1f, P0=%.3f)  Error: %f, Optimized α: %.7f, Optimized P: %.15f\n', ...
        initial_guesses_fractional(j, 1), initial_guesses_fractional(j, 2), ...
        errors_fractional(j), optimized_params_fractional(j, 1), optimized_params_fractional(j, 2));
   
    % Calculate efficiency gain
    if j <= length(errors_classical)
        efficiency_gain = abs(errors_classical(j) - errors_fractional(j)) / abs(errors_classical(j));
        fprintf('Efficiency Gain: %.6f          (%.2f%% improvement)\n', efficiency_gain, efficiency_gain * 100);
    end

end
