% Set the group number G.
G = 51;

% Define the range for n to cover all non-zero parts of x(n), h(n), and y(n).
% N_max is set to 3*G + 31 to ensure it includes all relevant indices.
N_max = 3*G + 31;
n = 0:N_max; % Create a vector n from 0 to N_max
L = length(n); % L is the total length of the vector n

% Define the input signal x(n): 1 for 0 <= n <= G + 10, 0 otherwise
x = zeros(1, L); % Initialize x as a vector of zeros with length L
x((n >= 0) & (n <= G + 10)) = 1; % Set x(n) = 1 for n in [0, G+10]

% Define the impulse response h(n) = x(n - G - 11): 1 for G + 11 <= n <= 2*G + 21, 0 otherwise
h = zeros(1, L); % Initialize h as a vector of zeros with length L
h((n >= G + 11) & (n <= 2*G + 21)) = 1; % Set h(n) = 1 for n in [G+11, 2*G+21]

% Compute the output y(n) using manual convolution.
% Since x(n) and h(n) are binary (0 or 1), y(n) is the number of k where x(k) = 1 and h(n - k) = 1.
y = zeros(1, L); % Initialize y as a vector of zeros with length L
for i = 1:L
    % For each i (corresponding to n = i - 1), find the range of k where x(k) = 1 and h(i - k) = 1
    % x(k) = 1 for 0 <= k <= G + 10, so k ranges from 1 to G + 11 in 1-based indexing
    % h(i - k) = 1 for G + 11 <= i - k <= 2*G + 21, so k <= i - (G + 11) and k >= i - (2*G + 21)
    low = max(1, i - (2*G + 21)); % Lower bound of k where h(i - k) could be 1
    high = min(G + 11, i - (G + 11)); % Upper bound of k where x(k) could be 1
    if high >= low
        y(i) = high - low + 1; % Number of overlapping k where both x(k) and h(i - k) are 1
    else
        y(i) = 0; % No overlap, so y(i) = 0
    end
end

% Plot the signals x(n), h(n), and y(n) using subplot and stem commands
figure; % Create a new figure for plotting
subplot(3,1,1); % First subplot: x(n)
stem(n, x); % Use stem to plot discrete-time signal x(n)
title('x(n)'); % Title of the plot
xlabel('n'); % Label for x-axis
ylabel('Amplitude'); % Label for y-axis

subplot(3,1,2); % Second subplot: h(n)
stem(n, h); % Use stem to plot discrete-time signal h(n)
title('h(n)'); % Title of the plot
xlabel('n'); % Label for x-axis
ylabel('Amplitude'); % Label for y-axis

subplot(3,1,3); % Third subplot: y(n)
stem(n, y); % Use stem to plot discrete-time signal y(n)
title('y(n)'); % Title of the plot
xlabel('n'); % Label for x-axis
ylabel('Amplitude'); % Label for y-axis