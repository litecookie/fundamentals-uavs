% Parameters
m = 0.420; % kg
g = 9.81; 
Cdf  = 0.97;
S = 0.01; % m^2
rho = 1.225; % kg/m^3
W = m*g; % N

% APC 8x6 SF Properties
% https://m-selig.ae.illinois.edu/props/volume-1/data/apcsf_8x6_static_2783rd.txt
dia = 0.2032; % m
radius = dia/2;
A = 4*pi*((radius)^2); %

num_cells = 3;
battery_life = 1500; % mA-hr
battery_voltage = 3.7; % V, assumed from propulsion lecture slide 53
Eb = num_cells * battery_voltage * (battery_life/1000) * 3600; % J

eta_m = 0.70;
eta_e = 0.80;

V = 0.5:0.5:20;

syms v;

solutions = zeros(size(V));
Pind = zeros(size(V));
Ptot = zeros(size(V));
T = zeros(size(V));
PtotbyV = zeros(size(V));
Df = zeros(size(V));

for i = 1:length(V)

    Df = (0.5 .* rho .* S .* Cdf .* (V.^2));

    T(i) = sqrt((W^2) + (Df(i)^2));
    alpha_D = atan(Df(i)/W);
    
    eq = v^4 + ((2 * V(i) * abs(sin(alpha_D)))*v^3) + (V(i)^2 * v^2) - (T(i)^2/(2 * rho * A)^2);

    % Solve for v
    sol = double(solve(eq, v));
    solutions(i) = sol(1); % Assuming there is a unique real solution

    Pind(i) = T(i) * solutions(i);
    
    sin_alpha_D = abs(sin(alpha_D));
    Ptot(i) = T(i) * (solutions(i) + (V(i) * sin_alpha_D));

    PtotbyV(i) = Ptot(i) / V(i);
end

% Plot of T vs V 
figure;
plot(V, T, 'b-', 'LineWidth', 2);
title('Thrust (T) vs Velocity(V)');
xlabel('V(m/s)');
ylabel('T(N)');
grid on;

% Plot of Pind vs V 
figure;
plot(V, Pind, 'g-', 'LineWidth', 2);
title('Induced Power (Pind) vs Velocity(V)');
xlabel('V(m/s)');
ylabel('Pind(W)');
grid on;

% Plot of Ptot vs V
figure;
plot(V, Ptot, 'Color', [1, 0.647, 0], 'Marker', 'o', 'LineWidth', 2);
title('Total Power (Ptot) vs Velocity(V)');
xlabel('V(m/s)');
ylabel('Ptot(W)');
grid on;

% Finding the minimum value of Ptot
[min_Ptot, min_Ptot_index] = min(Ptot);
min_Ptot_V = V(min_Ptot_index);

fprintf('Minimum Ptot value: %f at V = %f\n', min_Ptot, min_Ptot_V);

% Marking the minimum value on the Ptot graph
hold on;
scatter(min_Ptot_V, min_Ptot, 'ro', 'filled');
text(min_Ptot_V, min_Ptot, ['   Min Ptot: ', num2str(min_Ptot)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
hold off;

% Plot of Ptot/V vs V 
figure;
plot(V, PtotbyV, 'oc-', 'LineWidth', 2);
hold on;
title('Ptot/V vs Velocity(V)');
xlabel('V(m/s)');
ylabel('Ptot/V');
grid on;

% Finding the minimum value of Ptot/V
[min_PtotbyV, min_PtotbyV_index] = min(PtotbyV);
min_PtotbyV_V = V(min_PtotbyV_index);

fprintf('Minimum Ptot/V value: %f at V = %f\n', min_PtotbyV, min_PtotbyV_V);

% Marking the minimum value on the graph
scatter(min_PtotbyV_V, min_PtotbyV, 'ro', 'filled');
text(min_PtotbyV_V, min_PtotbyV, ['   Min Ptot/V: ', num2str(min_PtotbyV)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

hold off;

P = 9.351159671155774; %W % Ptot minimum value obtained from graph 
PbyV = 1.121846070222613; % Ptot/V minimum value obtained from graph

tmax = (Eb * eta_m * eta_e) / (P); %s 
fprintf('Max Endurance (tmax): %.2f s\n', tmax);

Rmax = (Eb * eta_m * eta_e) / (PbyV); %m
fprintf('Max Range (Rmax): %.2f m\n', Rmax);