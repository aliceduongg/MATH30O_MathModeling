function dydt = disease_model(t, y, params)
    % Parameters
    beta = params(1);
    gamma = params(2);
    delta = params(3);
    eta = params(4);
    alpha1 = params(5);
    alpha2 = params(6);
    
    % Variables
    S = y(1);
    I = y(2);
    U = y(3);
    W = y(4);
    R = y(5);
    
    % Equations
    dydt = zeros(5,1);
    dydt(1) = -beta*S*(I+U); % dS/dt
    dydt(2) = beta*S*(I+U) - (gamma+delta)*I; % dI/dt
    dydt(3) = delta*I - (eta+alpha1)*U; % dU/dt
    dydt(4) = gamma*I - (eta+alpha2)*W; % dW/dt
    dydt(5) = eta*W + eta*U; % dR/dt
end


% Parameters
beta = 4.44e-8;
gamma = 0.1142;
delta = 0.0285;
alpha1 = 1.5e-4;
alpha2 = 1.7826e-5;

% Different values of delta
eta_values = [0.1428, 0.05, 0.02];

% Time span
tspan = [0 300];

% Plotting I for different eta values
figure;
hold on;
for i = 1:length(eta_values)
    % Define delta
    eta = eta_values(i);
    
    % Parameters
    params = [beta, gamma, delta, eta, alpha1, alpha2];
    
    % Solve ODE
    [t, y] = ode45(@(t,y) disease_model(t, y, params), tspan, y0);
    
    % Plot I vs. time
    plot(t, y(:,3), 'LineWidth', 1.5);
end
hold off;
xlabel('Time');
ylabel('Asymptomatic Infected (U)');
legend('\eta = 0.1428', '\eta = 0.05', '\eta = 0.02');
title('Asymptomatic Infected (I) vs. Time for Different \eta Values');
grid on;

% Plotting U for different eta values
figure;
hold on;
for i = 1:length(eta_values)
    % Define delta
    eta = eta_values(i);

    % Parameters
    params = [beta, gamma, delta, eta, alpha1, alpha2];
    
    % Solve ODE
    [t, y] = ode45(@(t,y) disease_model(t, y, params), tspan, y0);
    
    % Plot U vs. time
    plot(t, y(:,5), 'LineWidth', 1.5);
end
hold off;
xlabel('Time');
ylabel('Recovered individuals (R)');
legend('\eta = 0.1428', '\eta = 0.05', '\eta = 0.02');
title('Unreported Symptomatic Infected (R) vs. Time for Different \eta Values');
grid on;

% Plotting W for different eta values
figure;
hold on;
for i = 1:length(eta_values)
    % Define delta
    eta = eta_values(i);
    
    % Parameters
    params = [beta, gamma, delta, eta, alpha1, alpha2];
    
    % Solve ODE
    [t, y] = ode45(@(t,y) disease_model(t, y, params), tspan, y0);
    
    % Plot W vs. time
    plot(t, y(:,4), 'LineWidth', 1.5);
end
hold off;
xlabel('Time');
ylabel('Reported Symptomatic Infected (W)');
legend('\eta = 0.1428', '\eta = 0.05', '\eta = 0.02');
title('Reported Symptomatic Infected (W) vs. Time for Different \eta Values');
grid on;
