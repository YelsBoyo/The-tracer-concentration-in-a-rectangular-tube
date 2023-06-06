% Parameters
L = 110000; % Length of the channel (m)
N = 100; % Number of grid points
dx = L / N; % Grid spacing (m)
dt = 10; % Time step (s)
T = 36000; % Total simulation time (s)

% Constants
H = 0.4; % Channel height (m)
B = 10; % Channel width (m)
S0 = 0.0004; % Initial concentration
v = 1; % Velocity (m/s)
D = 0.1; % Dispersion coefficient (m^2/s)

% Calculate the Courant number
C = v * dt / dx;

% Initialize concentration arrays
C_uw = zeros(N, T / dt);
C_lw = zeros(N, T / dt);
C_quickest = zeros(N, T / dt);

% Set initial concentration
C_uw(:, 1) = S0;
C_lw(:, 1) = S0;
C_quickest(:, 1) = S0;

% Numerical simulation
for t = 2:T/dt
    % Upwind method
    for i = 2:N-1
        C_uw(i, t) = C_uw(i, t-1) - C * (C_uw(i, t-1) - C_uw(i-1, t-1));
    end
    
    % Lax-Wendroff method
    for i = 2:N-1
        C_lw(i, t) = C_lw(i, t-1) - C/2 * (C_lw(i+1, t-1) - C_lw(i-1, t-1))...
            + C^2/2 * (C_lw(i+1, t-1) - 2 * C_lw(i, t-1) + C_lw(i-1, t-1));
    end
    
    % QUICKEST method
    for i = 5:N-4 % Adjusted loop range
        C_quickest(i, t) = C_quickest(i, t-1) - C * (9 * C_quickest(i, t-1) - 116 * C_quickest(i-1, t-1)...
            + 213 * C_quickest(i-2, t-1) - 48 * C_quickest(i-3, t-1) + 8 * C_quickest(i-4, t-1)) / 160;
    end
    
    % Boundary conditions
    C_uw(1, t) = C_uw(2, t);
    C_uw(N, t) = C_uw(N-1, t);
    
    C_lw(1, t) = C_lw(2, t);
    C_lw(N, t) = C_lw(N-1, t);
    
    C_quickest(1, t) = C_quickest(2, t);
    C_quickest(2, t) = C_quickest(3, t);
    C_quickest(N, t) = C_quickest(N-1, t);
    C_quickest(N-1, t) = C_quickest(N-2, t);
end

% Visualization
x = dx * (0:N-1);

figure;
hold on;
plot(x, C_uw(:, end), 'r-', 'LineWidth', 1.5);
plot(x, C_lw(:, end), 'g--', 'LineWidth', 1.5);
plot(x, C_quickest(:, end), 'b-.', 'LineWidth', 1.5);
hold off;
xlabel('Distance (m)');
ylabel('Concentration');
legend('Upwind', 'Lax-Wendroff', 'QUICKEST');
title('Concentration at t = T');

figure;
mesh(1:T/dt, x, C_uw);
xlabel('Time step');
ylabel('Distance (m)');
zlabel('Concentration');
title('Concentration using Upwind method');

figure;
mesh(1:T/dt, x, C_lw);
xlabel('Time step');
ylabel('Distance (m)');
zlabel('Concentration');
title('Concentration using Lax-Wendroff method');

figure;
mesh(1:T/dt, x, C_quickest);
xlabel('Time step');
ylabel('Distance (m)');
zlabel('Concentration');
title('Concentration using QUICKEST method');
