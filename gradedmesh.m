% Solves a fractional diffusion equation with a specified fractional order alpha
% Inputs:
%   alpha: Fractional order of the time derivative (0 < alpha < 1)
%   T: Final time
%   X_val: Spatial domain length (0 to X_val)
%   M: Number of spatial intervals
%   N: Number of time intervals
%   r: Parameter controlling time mesh grading
% Outputs:
%   U: Solution matrix (N+1 x M+1), where U(n,m) approximates u(t_n, x_m)
%   t: Time mesh (N+1 x 1)
%   x: Spatial mesh (M+1 x 1)
function [U, t, x] = solve_fractional_diffusion(alpha, T, X_val, M, N, r)
    % Generate time mesh using a graded mesh formula t_k = T * (k/N)^r
    t = zeros(N+1,1);
    for k=0:N
        t(k+1) = T * (k/N)^r; % Non-uniform time steps, concentrated near t=0 for r>1
    end
    
    % Generate spatial mesh with uniform spacing
    x = linspace(0, X_val, M+1); % x_m = m * h, where h = X_val / M
    h = x(2) - x(1); % Spatial step size
    
    % Initialize solution matrix U
    U = zeros(N+1, M+1); % U(n+1,m+1) stores u(t_n, x_m)
    U(1,:) = sin(x); % Initial condition: u(x,0) = sin(x)
    
    % Precompute Gamma(2 - alpha) for fractional derivative coefficients
    Gamma_2_minus_alpha = gamma(2 - alpha);
    
    % Construct spatial operator matrix A for interior points (m=2 to M)
    % A represents the discretized second derivative plus a reaction term
    A = zeros(M-1, M-1); % Tridiagonal matrix for interior points
    for i=1:M-1
        A(i,i) = 2/h^2 + (1 + x(i+1)); % Diagonal: 2/h^2 + (1 + x_m)
        if i > 1
            A(i,i-1) = -1/h^2; % Sub-diagonal: -1/h^2
        end
        if i < M-1
            A(i,i+1) = -1/h^2; % Super-diagonal: -1/h^2
        end
    end
    
    % Time stepping loop to compute solution at each time level
    for n=1:N
        j = n+1; % Index for current time step (t_j = t(n+1))
        tau_n = t(j) - t(j-1); % Current time step size
        
        % Compute coefficients b_k for the fractional derivative approximation
        b = zeros(n,1); % b_k for k=0 to n-1
        for k=0:n-1
            tau_kp1 = t(k+2) - t(k+1); % Time step size at t_(k+1)
            if tau_kp1 == 0
                continue; % Skip if time step is zero (avoid division by zero)
            end
            % Compute terms for the fractional derivative
            term1 = (t(j) - t(k+1))^(1-alpha);
            term2 = (t(j) - t(k+2))^(1-alpha);
            b(k+1) = 1 / Gamma_2_minus_alpha / tau_kp1 * (term1 - term2);
        end
        
        % Compute coefficients a_vec for past solution terms
        a_vec = zeros(n+1,1); % Coefficients for U(1,m) to U(n+1,m)
        a_vec(1) = -b(1); % Coefficient for U(1,m) (initial condition)
        for l=2:n+1
            if l <= n
                a_vec(l) = b(l-1) - b(l); % Coefficients for U(2,m) to U(n,m)
            else
                a_vec(n+1) = b(n); % Coefficient for U(n+1,m)
            end
        end
        
        % Compute historical contribution P_j_int using past solutions
        U_hist = U(1:n,:); % Solutions from t_0 to t_(n-1)
        P_j_full = a_vec(1:n)' * U_hist(:,2:M); % Sum of weighted past solutions at interior points
        
        % Compute source term f_j_int at interior points
        f_j_int = (x(2:M).*(pi - x(2:M)).*(1 + t(j)^4) + t(j)^2)'; % f(x,t) = x*(pi-x)*(1+t^4) + t^2
        
        % Coefficient of U(n+1,m) in the fractional derivative
        c_j = a_vec(n+1); % Coefficient for the current solution
        
        % Solve the linear system for interior points
        matrix = c_j * eye(M-1) + A; % System matrix: c_j * I + A
        rhs = f_j_int - P_j_full'; % Right-hand side: f_j - P_j
        U_j_int = matrix \ rhs; % Solve for U(n+1,2:M)
        
        % Update solution at current time step
        U(j,2:M) = U_j_int'; % Assign interior points
        U(j,1) = 0; % Boundary condition: u(0,t) = 0
        U(j,M+1) = 0; % Boundary condition: u(X_val,t) = 0
    end
end

% Computes the maximum difference between coarse and fine mesh solutions
% Inputs:
%   U_coarse, U_fine: Solution matrices on coarse and fine meshes
%   t_coarse, t_fine: Time meshes for coarse and fine grids
%   x_coarse, x_fine: Spatial meshes for coarse and fine grids
% Output:
%   D: Maximum absolute difference between coarse and fine solutions
function D = compute_two_mesh_difference(U_coarse, U_fine, t_coarse, t_fine, x_coarse, x_fine)
    [Nc, Mc] = size(U_coarse); % Dimensions of coarse solution
    [Nf, Mf] = size(U_fine); % Dimensions of fine solution
    
    % Validate that fine mesh is twice as fine as coarse mesh
    if Mf ~= 2*Mc - 1 || Nf ~= 2*Nc - 1
        error('Fine mesh must be twice the coarse mesh in both dimensions');
    end
    
    max_diff = 0; % Initialize maximum difference
    for i=1:Nc
        for j=1:Mc
            % Get coarse grid point coordinates
            t_c = t_coarse(i);
            x_c = x_coarse(j);
            
            % Find nearest corresponding point in fine grid
            t_f_idx = find(abs(t_fine - t_c) == min(abs(t_fine - t_c)), 1);
            x_f_idx = find(abs(x_fine - x_c) == min(abs(x_fine - x_c)), 1);
            
            % Compute absolute difference at this point
            diff_val = abs(U_coarse(i,j) - U_fine(t_f_idx, x_f_idx));
            if diff_val > max_diff
                max_diff = diff_val; % Update maximum difference
            end
        end
    end
    D = max_diff; % Return maximum difference
end

% Main script for generating plot and convergence table
clear all;
alpha = 0.6; % Fractional order
T = 1; % Final time
X_val = pi; % Spatial domain [0, pi]
M_plot = 100; % Spatial intervals for plot
N_plot = 100; % Time intervals for plot
r_plot = (2 - alpha)/alpha; % Time mesh grading parameter for plot

% Plot the solution as a surface
[U_plot, t_plot, x_plot] = solve_fractional_diffusion(alpha, T, X_val, M_plot, N_plot, r_plot);
figure(1);
[ X_grid, Y_grid ] = meshgrid(x_plot, t_plot'); % Create 2D grid for plotting
surf(X_grid, Y_grid, U_plot'); % Plot surface
xlabel('x'); % Label x-axis
ylabel('t'); % Label y-axis
zlabel('u(x,t)'); % Label z-axis
title('Solution u(x,t) for \alpha = 0.6'); % Title

% Generate convergence table
r_values = [1, (2 - alpha)/(2*alpha), (2 - alpha)/alpha, 2*(2 - alpha)/alpha]; % Time mesh grading parameters
M_list = [64; 128; 256; 512; 1024; 2048]; % Mesh sizes
%M_list = [4]; % Uncomment for testing with a small mesh
fprintf('Table for alpha = 0.6:\n');
fprintf('r\tN=M=64\tN=M=128\tN=M=256\tN=M=512\tN=M=1024\n');

% Loop over each grading parameter
for k=1:length(r_values)
    r = r_values(k);
    fprintf('r=%g\n', r);
    D_list = zeros(1, length(M_list)-1); % Store differences between consecutive meshes
    
    % Compute solutions on consecutive meshes and their differences
    for j=1:length(M_list)-1
        M_coarse = M_list(j); % Coarse mesh size
        N_coarse = M_coarse; % Same number of time and spatial intervals
        M_fine = M_list(j+1); % Fine mesh size
        N_fine = M_fine;
        
        % Solve on coarse and fine meshes
        [U_coarse, t_coarse, x_coarse] = solve_fractional_diffusion(alpha, T, X_val, M_coarse, N_coarse, r);
        [U_fine, t_fine, x_fine] = solve_fractional_diffusion(alpha, T, X_val, M_fine, N_fine, r);
        
        % Compute maximum difference
        D_list(j) = compute_two_mesh_difference(U_coarse, U_fine, t_coarse, t_fine, x_coarse, x_fine);
    end
    
    % Compute convergence orders
    p_list = zeros(1, length(M_list)-2); % Store convergence orders
    for i=1:length(M_list)-2
        p_list(i) = log2(D_list(i) / D_list(i+1)); % Order p = log2(D_coarse / D_fine)
    end
    
    % Print results
    fprintf('D: ');
    for i=1:length(D_list)
        fprintf('%e\t', D_list(i)); % Print differences
    end
    fprintf('\n');
    fprintf('p: ');
    for i=1:length(p_list)
        fprintf('%f\t', p_list(i)); % Print convergence orders
    end
    fprintf('\n\n');
end