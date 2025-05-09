% MATLAB code to solve fractional differential equation using RBFs 
% Equation: D^0.8 y(x) = (Gamma(3)/Gamma(2.2)) * x^1.2, y(0) = 0
% Exact solution: y(x) = x^2

% Define RBF function with stability improvements
function phi = rbf(type, r, beta, nu)
    switch type
        case 'gaussian'
            phi = exp(-r^2 / 2);
        case 'multiquadric'
            phi = (1 + r^2)^(beta/2);
        case 'matern_sobolev'
            c = 10; % Scaling factor for stability
            phi = (c*r)^nu * besselk(nu, c*r); 
        case 'powers'
            phi = r^beta;
        case 'thin_plate'
            phi = r^(2*beta) * log(r + 1e-20);
        otherwise
            error('Unknown RBF type');
    end
end

% Modified phi_j(x) with stability checks
function p = phi_j(x, j, type, beta, nu, x_tilde)
    r = abs(x - x_tilde(j));
    % Special handling for Matérn-Sobolev at r=0
    if strcmp(type, 'matern_sobolev') && r < 1e-20
        r = 1e-20; % Prevent singularity
    end
    p = rbf(type, r, beta, nu);
end

% Grünwald-Letnikov fractional derivative 
function df = gl_fd(f, alpha, x, h)
    m = floor(x / h);
    df = 0;
    for k = 0:m
        binom = 1;
        for j = 1:k
            binom = binom * (alpha - j + 1) / j;
        end
        t = x - k * h;
        if t < 0
            continue
        end
        df = df + (-1)^k * binom * f(t);
    end
    df = df / h^alpha;
end

% Main code with stability improvements
n = 20; % Number of collocation points
x = linspace(0,50,n);
x_tilde = x;
alpha = 0.8; % Fractional order
h = 1e-4; % Step size
reg_param = 1e-12; % Increased regularization for stability

% Updated RBF parameters for Matérn-Sobolev
types = {'gaussian', 'multiquadric', 'matern_sobolev', 'powers', 'thin_plate'};
betas = [1, 1, 1.5, 3, 2];
nus = [0, 0, 1.5, 0, 0]; 

% Compute all approximated solutions
y_approx_all = cell(1,5);
for type_idx = 1:5
    type = types{type_idx};
    beta_val = betas(type_idx);
    nu_val = nus(type_idx);
    
    % Initialize matrix C
    C = zeros(n,n);
    
    % Boundary condition row
    for j = 1:n
        C(1,j) = phi_j(x(1), j, type, beta_val, nu_val, x_tilde);
    end
    
    % Fractional derivative rows
    for i = 2:n
        for j = 1:n
            phi_j_func = @(t) phi_j(t, j, type, beta_val, nu_val, x_tilde);
            C(i,j) = gl_fd(phi_j_func, alpha, x(i), h);
        end
    end
    
    % Regularization and solve
    C = C + reg_param * eye(n);
    G = zeros(n, 1);
    
    %G(2:n) = (gamma(3)/gamma(2.2)) * x(2:n).^1.2;
    G(2:n) = x(2:n) .* exp(-x(2:n));

    try
        lambda = C \ G;
    catch ME
        warning(['Failed to solve for ' type ' RBF: ' ME.message]);
        y_approx_all{type_idx} = NaN(size(linspace(0,1,200)));
        continue;
    end
    
    % Evaluate solution
    x_plot = linspace(0,1,200);
    y_approx = zeros(size(x_plot));
    for k = 1:length(x_plot)
        for j = 1:n
            y_approx(k) = y_approx(k) + lambda(j) * phi_j(x_plot(k), j, type, beta_val, nu_val, x_tilde);
        end
    end
    y_approx_all{type_idx} = y_approx;
    
    % Error calculation
    y_exact = x_plot.^2;
    if all(~isnan(y_approx))
        error_inf = max(abs(y_exact - y_approx));
        fprintf('L-infinity error for %-15s: %.6e\n', type, error_inf);
    end
end

% Plot results (unchanged)
figure;
x_plot = linspace(0,1,200);
y_exact = x_plot.^2;
plot(x_plot, y_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Exact');
hold on;
colors = ['r', 'g', 'c', 'm', 'k'];
line_styles = ['--', '-.', ':', '-', '--'];
for type_idx = 1:5
    if ~isnan(y_approx_all{type_idx}(1))
        type = types{type_idx};
        plot(x_plot, y_approx_all{type_idx}, [colors(type_idx) line_styles(type_idx)],...
            'LineWidth', 2, 'DisplayName', [type ' Approximation']);
    end
end
hold off;
legend('Location', 'northwest');
title('Solution using different RBFs for D^{0.8} y(x) = (Γ(3)/Γ(2.2)) x^{1.2}');
xlabel('x');
ylabel('y(x)');
grid on;
