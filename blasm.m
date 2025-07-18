
% MATLAB script to solve the Blasius equation using shooting method with RK4
% The Blasius equation is f''' + (1/2) * f * f'' = 0
% Boundary conditions: f(0) = 0, f'(0) = 0, f'(eta -> inf) = 1

% Define parameters
eta_max = 10; % Approximate infinity
n = 1000; % Number of points
deta = eta_max / (n - 1); % Step size
eta = 0:deta:eta_max; % Eta grid

% Initial guess for f''(0) (to be adjusted via shooting)
f2_initial = 0.332; % Initial guess for f''(0)

% Tolerance for convergence of f'(eta_max) to 1
tol = 1e-6;
max_iter = 100;

% Shooting method
f2 = f2_initial;
df2 = 0.01; % Perturbation for secant method
for iter = 1:max_iter
    % Solve ODE with current f2
    [f, fp, fpp] = solve_blasius(eta, f2);
    
    % Check boundary condition f'(eta_max) = 1
    error = fp(end) - 1;
    if abs(error) < tol
        fprintf('Converged after %d iterations, f''''(0) = %.6f\n', iter, f2);
        break;
    end
    
    % Perturb f2 and solve again
    [~, fp_perturb, ~] = solve_blasius(eta, f2 + df2);
    error_perturb = fp_perturb(end) - 1;
    
    % Secant method to update f2
    f2_new = f2 - error * (df2 / (error_perturb - error));
    f2 = f2_new;
end

% Plot results
figure;
subplot(2,1,1);
plot(eta, f, 'b-', eta, fp, 'r--', eta, fpp, 'g-.');
legend('f (displacement)', 'f'' (velocity)', 'f'''' (shear)', 'Location', 'best');
xlabel('\eta'); ylabel('f, f'', f''''');
title('Blasius Solution');
grid on;

subplot(2,1,2);
plot(eta, fp, 'r-');
xlabel('\eta'); ylabel('f'' (velocity)');
title('Velocity Profile');
grid on;

% Function to solve the ODE system using RK4
function [f, fp, fpp] = solve_blasius(eta, f2)
    n = length(eta);
    deta = eta(2) - eta(1);
    
    % Initialize arrays
    f = zeros(1, n);
    fp = zeros(1, n);
    fpp = zeros(1, n);
    
    % Initial conditions
    f(1) = 0; % f(0) = 0
    fp(1) = 0; % f'(0) = 0
    fpp(1) = f2; % f''(0) = initial guess
    
    % RK4 integration
    for i = 1:n-1
        % Current state
        y = [f(i); fp(i); fpp(i)];
        
        % RK4 steps
        k1 = blasius_ode(y);
        k2 = blasius_ode(y + 0.5*deta*k1);
        k3 = blasius_ode(y + 0.5*deta*k2);
        k4 = blasius_ode(y + deta*k3);
        
        % Update state
        y = y + (deta/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        % Store results
        f(i+1) = y(1);
        fp(i+1) = y(2);
        fpp(i+1) = y(3);
    end
end

% Function defining the ODE system
function dydeta = blasius_ode(y)
    f = y(1);
    fp = y(2);
    fpp = y(3);
    dydeta = [fp; fpp; -0.5 * f * fpp];
end
