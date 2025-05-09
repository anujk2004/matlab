% RK 2nd Order
clc; clear; close all;

%Start time
tic;

% DE to put below
f = @(x, y)(1/(x+1)*y^2 ) ; 

% Define Symbolic Variables
syms y(x)

% Defining the Differential Equation
ode = diff(y, x) == f(x,y);

% Define Initial Condition
y0 = 1;
cond = y(0) == y0;

% to solve the Differential Equation 
sol = dsolve(ode, cond);

% Convert to MATLAB Function for Numeric Evaluation
actual_solution = matlabFunction(sol);


% Input parameters
a = 0;    % Start
b = 1;    % End 
N = 10;   % Number of subintervals

h = (b - a) / N;

x_val_array = a:h:b;

y_approx_array = zeros(size(x_val_array));
y_approx_array(1) = 1; % Initial condition


% The algo is below
for i = 1:N
    K1 = h * f(x_val_array(i), y_approx_array(i));
    K2 = h * f(x_val_array(i) + h, y_approx_array(i) + K1);

    y_approx_array(i+1) = y_approx_array(i) + (K1+K2)/2;
     
end


% Compute actual solution at x points
y_actual_array = actual_solution(x_val_array);

% Display results
disp('X values:');
disp(x_val_array);
disp('RK method 2:');
disp(y_approx_array);
disp('Actual Solution:');
disp(y_actual_array);

% Stop timing
elapsed_time = toc;  
fprintf('Execution Time: %.6f seconds\n', elapsed_time);

% Plot
figure;
plot(x_val_array, y_approx_array, '-o', 'LineWidth', 1.5, 'DisplayName', 'RK method of 2nd Order');
hold on;
plot(x_val_array, y_actual_array, '-x', 'LineWidth', 1.5, 'DisplayName', 'Actual Solution');
xlabel('x');
ylabel('y');
legend;
grid on;
title('RK method 2''s Method vs. Actual Solution');

% Compute Errors for Order of Convergence
N_values = [N, 2*N, 4*N, 8*N,16*N];  % Different N values for convergence analysis


errors = zeros(size(N_values));

for j = 1:length(N_values)
    N = N_values(j);
    h = (b - a) / N;   % Step size
    x_val_array = a:h:b;
    
    % Initialize y_approx array
    y_approx_array = zeros(size(x_val_array));
    y_approx_array(1) = y0; % Initial condition

    % Euler's Method
    for i = 1:N
         K1 = h * f(x_val_array(i), y_approx_array(i));
         K2 = h * f(x_val_array(i) + h, y_approx_array(i) + K1);

         y_approx_array(i+1) = y_approx_array(i) + (K1+K2)/2;
    end

    % Compute actual solution
    y_actual_array = actual_solution(x_val_array);

    % Compute error (norm of the difference)
    errors(j) = max(abs(y_actual_array - y_approx_array));
    
    % Display results
    fprintf('For N = %d, Max Error = %.6e\n', N, errors(j));
end

% Compute Order of Convergence
fprintf('\n Order of Convergence (p):\n');
for j = 2:length(N_values)
    h1 = (b - a) / N_values(j-1);  
    h2 = (b - a) / N_values(j);  
    p = log(errors(j-1) / errors(j)) / log(h1/h2);
    fprintf('Between N=%d (h=%.6f)and N=%d (h=%.6f), p = %.15f\n', N_values(j-1), h1, N_values(j), h2, p);
end

% Stop timing
elapsed_time = toc; 
fprintf( '\nExecution Time: %.6f seconds\n',elapsed_time );



