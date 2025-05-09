% main.m
% Check if exact solution exists and run the appropriate solver
clear all; clc;




% Define exact solution (comment out if not available)
u_exact = @(x, t) t.^2 .* sin(pi*x);  % Example exact solution

% --- Check if Exact Solution Exists ---
if exist('u_exact', 'var') && isa(u_exact, 'function_handle')
    fprintf('Exact solution detected. Running solgiven.m...\n');
    solgiven;  % Run solgiven.m (computes errors using exact solution)
else
    fprintf('No exact solution detected. Running gradedmesh.m...\n');
    gradedmesh; % Run gradedmesh.m (uses two-mesh difference)
end