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

% Initial conditions
S0 = 11.081e6;
I0 = 3.62;
U0 = 0.2;
W0 = 4.13;
R0 = 0;

% Parameters
beta = 4.44e-8;
gamma = 0.1142;
delta = 0.0285;
eta = 1/7;
alpha1 = 1.5e-4;
alpha2 = 1.7826e-5;

params = [beta, gamma, delta, eta, alpha1, alpha2];

% Combine initial conditions
y0 = [S0; I0; U0; W0; R0];

% Time span
tspan = [0 100]; % Example time span, adjust as needed

% Solve ODE
[t, y] = ode45(@(t,y) disease_model(t, y, params), tspan, y0);


% Plotting S vs. I
figure;
plot(y(:,2), y(:,1), 'b');
xlabel('Asymptomatic Infected (I)');
ylabel('Susceptible (S)');
title('Susceptible vs. Asymptomatic Infected');
grid on;

% Plotting U vs. I
figure;
plot(y(:,2), y(:,3), 'g');
xlabel('Asymptomatic Infected (I)');
ylabel('Unreported Symptomatic Infected (U)');
title('Unreported Symptomatic vs. Asymptomatic Infected');
grid on;

% Plotting W vs. I
figure;
plot(y(:,2), y(:,4), 'm');
xlabel('Asymptomatic Infected (I)');
ylabel('Reported Symptomatic Infected (W)');
title('Reported Symptomatic vs. Asymptomatic Infected');
grid on;

% Plotting R vs. I
figure;
plot(y(:,2), y(:,5), 'k');
xlabel('Asymptomatic Infected (I)');
ylabel('Recovered (R)');
title('Recovered vs. Asymptomatic Infected');
grid on;