SSSprop = 0.2027; %m^2
S = 0.55; %m^2
r = sqrt(Sprop/pi); %m
d = 2 * r; %m
T = 13.1 * 9.81; %N
W = T; %steady level cruise flight
rho = 1.225; %kg/m^3
rpm = 12500; 
n = 12500/60; %rps
Ct = (T / (rho * (n^2) * (d^4)));
Jme = 0.6920;
V = Jme * n * d; %m^2/s
Cq = 0.0050;
Cp = 2 * pi * Cq;
P = (Cp * rho * (n^3) * (d^5)); %W
sfc = ((0.6651/(1000*3600)) * 9.81); %N/s/W 


eta = ((Ct * Jme) / Cp); 
eta_check = ((T * V) / P);

e = 0.9;
b = 2.8956; %m
c = 0.18994; %m
AR = ((b^2)/S);
K = (1 / (pi * e * AR));
cdo = 0.03;
W0 = 13.1 * 9.81; %N
W1 = 9.1 * 9.81; %N

%Range and Endurance expressed in km and hour
Rmax = ((eta / sfc) * (1 / (2 * (sqrt(K * cdo)))) * log(W0/W1));
Emax =  (eta / sfc) * (sqrt(2 * rho * Sprop)) * ((1 / (4 * cdo))*(((3 * cdo) / K)^0.75)) * ((1/(sqrt(W1))) - (1/(sqrt(W0))));
Emax_hour = (Emax/3600);
Rmax_km = (Rmax/(10^3));

V_given = 19.44;
%Angle of Attack for Cruise Flight
CLalpha = 3.45;

alphacruise = ((((2*W)/(rho*(V_given^2)*S))) / CLalpha);

% Max Climb Rate 
Vgammamax = (4*K*(W/S))/(rho*((eta*P)/W));
singamma = ((eta*P)/W) - ((0.5 * rho * (Vgammamax^2)) * ((cdo)/(W/S))) - ((2*K*(W/S))/(rho*(V_given^2)));
Clprmin = sqrt((3*cdo)/K);
Vprmin = sqrt((2*W)/(rho*S*Clprmin));
RbyCmax = singamma * Vprmin;

% Rate of Climb - Theoretical - Ascend to 1100m.
RbyC = ((eta*P)/W) - ((0.5 * rho * (V_given^3)) * ((cdo)/(W/S))) - ((2*K*(W/S))/(rho*V_given));

%Yaw Rate for turn over radius of 250m - In Radians
yawrate = ((V_given)/(250));

% For Q7, edit each of u, v, w in the Va equation in aircraftdynamics.m -
% Based on simulation observations, write down theoretical answers

% Simulation Calculations
% airspeed velocity value in x-direction: 

Vx_simulation = 19.0635;
% Alpha Cruise
alphacruisesimulation = ((((2*W)/(rho*(Vx_simulation^2)*S))) / CLalpha);

% Rate of Climb
RbyCsimulation = ((eta*P)/W) - ((0.5 * rho * (Vx_simulation^3)) * ((cdo)/(W/S))) - ((2*K*(W/S))/(rho*Vx_simulation));

%Yaw Rate for turn over radius of 250m - In Radians
yawratesimulation = ((Vx_simulation)/(250));

error1 = round((abs(alphacruise - alphacruisesimulation)*100/(alphacruisesimulation)), 2);
error2 = round((abs(yawrate - yawratesimulation)*100/(alphacruisesimulation)), 2);
error3 = round((abs(RbyC - RbyCsimulation)*100/(alphacruisesimulation)), 2);

fprintf("\n\t\t\t Theoretical Value   Simulation Value    Rate of Error\n")
fprintf("Alpha Cruise %d \t     %d \t    %d\n", alphacruise, alphacruisesimulation, error1);
fprintf("Yaw Rate     %d \t     %d \t    %d\n", yawrate, yawratesimulation, error2);
fprintf("Alpha Cruise %d \t     %d \t    %d\n", RbyC, RbyCsimulation, error3);
