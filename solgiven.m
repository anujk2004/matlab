% Solves a time-fractional PDE at a fixed spatial point using the L1 scheme
% Computes numerical solutions for different fractional orders (alpha) and mesh sizes
% Outputs a table of maximum nodal errors and convergence rates at x = π/2

% Define parameters
alpha_values = [0.2, 0.4, 0.6, 0.8]; % Fractional orders of the time derivative
N_values = 4 * 2.^(0:6); % Temporal mesh sizes: [4, 8, 16, 32, 64, 128, 256]
T = 1; % Final time
x = pi/2; % Fixed spatial point for evaluation (x = π/2)
a = pi^2; % Coefficient of the solution term in the PDE (u term)

% Preallocate results array to store errors and convergence rates
% Dimensions: (number of alpha values) x (number of N values) x (2 for error and rate)
results = zeros(length(alpha_values), length(N_values), 2);

% Loop over each fractional order alpha
for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx); % Current fractional order
    errors = zeros(size(N_values)); % Array to store errors for each N
    
    % Loop over each mesh size
    for idx = 1:length(N_values)
        N = N_values(idx); % Number of temporal intervals
        M = N; % Number of spatial intervals (set equal to temporal, though not used directly)
        h = pi / M; % Spatial step size (domain: x ∈ [0, π], not used since x is fixed)
        tau = T / N; % Temporal step size
        
        % Define exact solution and source term
        % Exact solution: u(x,t) = t^2 * sin(πx)
        u_exact = @(x, t) t.^2 .* sin(pi*x);
        % Source term: f(x,t) = (2/Γ(3-α)) * t^(2-α) * sin(πx) + π^2 * t^2 * sin(πx)
        f = @(x, t) (2/gamma(3 - alpha)) * t.^(2 - alpha) .* sin(pi*x) + pi^2 * t.^2 .* sin(pi*x);
        
        % Compute L1 scheme coefficients for fractional derivative approximation
        b = zeros(N, 1); % Coefficients b_l for l=1 to N
        for l = 1:N
            b(l) = l^(1 - alpha) - (l-1)^(1 - alpha); % b_l = l^(1-α) - (l-1)^(1-α)
        end
        % Compute differences of b_l coefficients for historical sum
        if N > 1
            diff_b = b(2:N) - b(1:N-1); % diff_b(l) = b(l+1) - b(l) for l=1 to N-1
        else
            diff_b = []; % Empty for N=1, as no historical terms exist
        end
        
        % Initialize solution array for u at fixed x
        U = zeros(N, 1); % U(n) approximates u(x, t_n) at x = π/2
        
        % Time-stepping loop to compute solution at each time level
        for n = 1:N
            t_n = n * tau; % Current time t_n = n * τ
            f_n = f(x, t_n); % Source term at (x, t_n)
            
            % Compute historical sum for L1 scheme (fractional derivative)
            if n == 1
                sum_term = 0; % No historical terms at first time step
            else
                % Compute sum of past terms: sum_{l=1}^{n-1} (b(n-l+1) - b(n-l)) * U(l)
                sum_indices = n-1:-1:1; % Reverse indices for correct pairing
                sum_term = dot(U(1:n-1), diff_b(sum_indices));
            end
            
            % Compute L1 scheme coefficient for current time step
            coeff_L1 = b(1) / (gamma(2 - alpha) * tau^alpha); % Coefficient for U(n)
            
            % Generalized denominator including coefficient a
            denominator = coeff_L1 + a;
            
            % Update solution using L1 scheme
            % U(n) = (f_n - historical_sum / (Γ(2-α) * τ^α)) / denominator
            U(n) = (f_n - sum_term / (gamma(2 - alpha) * tau^alpha)) / denominator;
        end
        
        % Compute error at all time points for fixed x
        t = tau * (1:N)'; % Time grid t_n = n * τ
        exact_values = u_exact(x, t); % Exact solution at x = π/2, t = t_n
        errors(idx) = max(abs(U - exact_values)); % Maximum nodal error
    end
    
    % Calculate convergence rates
    rates = zeros(length(N_values)-1, 1); % Convergence rates for consecutive meshes
    for idx = 1:length(N_values)-1
        rates(idx) = log2(errors(idx)/errors(idx+1)); % Rate = log2(E_coarse / E_fine)
    end
    
    % Store errors and rates in results array
    results(a_idx, :, 1) = errors; % Store errors
    results(a_idx, 1:end-1, 2) = rates; % Store rates (no rate for last N)
end

% Display results in a formatted table
fprintf('\nResults for Example 3.3 (Maximum Nodal Errors and Rates)\n');
fprintf('Alpha   N       Error         Rate\n');
for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx); % Current fractional order
    for idx = 1:length(N_values)
        error_val = results(a_idx, idx, 1); % Error for current N
        if idx == 1
            rate_str = '  -  '; % No rate for first mesh
        else
            rate_val = results(a_idx, idx-1, 2); % Convergence rate
            rate_str = sprintf('%7.3f', rate_val);
        end
        % Print alpha, N, error, and rate
        fprintf('%.1f   %6d  %.3e  %s\n', alpha, N_values(idx), error_val, rate_str);
    end
    fprintf('\n'); % Blank line between alpha values
end