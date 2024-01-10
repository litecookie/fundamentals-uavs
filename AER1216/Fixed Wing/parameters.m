% AER1216 Fall 2023 
% Fixed Wing Project Code
%
% parameters.m
%
% Initialization file which generates and stores all required data into the 
% structure P, which is then stored in the workspace. Simulink model runs 
% this file at the start of every simulation. Code structure adapted from
% Small Unmanned Aircraft: Theory and Practice by R.W. Beard and T. W. 
% McLain. 
% 
% Inputs: 
% N/A
%
% Outputs:
% P                 structure that contains all aerodynamic, geometric, and
%                   initial condition data for the aircraft and simulation.
%
% Last updated: Pravin Wedage 2023-10-26

%% TA NOTE
% An easy way to store parameters for use in simulink is through the use of
% a structure. For example, P.g = 9.81 stores the value of gravitational
% acceleration in the field g that is contained within the structure P.
% Anytime P is called anywhere in the simulation code, the value of P.g is
% accessible. 

%% Parameter Computation
% Initial Conditions
clear all
% compute trim conditions            
P.Va0 = 20;         % initial airspeed (also used as trim airspeed)
P.Va_trim = 20; 
P.Va = P.Va_trim;


P.gravity = 9.81;
P.g = 9.81; 

% Aerosonde UAV Data
% physical parameters of airframe

% aerodynamic coefficients

% Control Input limits 
P.delta_e_max = deg2rad(45); % assumed symmetric
P.delta_a_max = deg2rad(45); 
P.delta_r_max = deg2rad(25);

% Initial Conditions % connects with aircraft_dynamics.m, do not modify
% structure field names
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -1000;  % initial Down position (negative altitude)
P.u0     = P.Va0; % initial velocity along body x-axis
P.v0     = 0;  % initial velocity along body y-axis
P.w0     = 0;  % initial velocity along body z-axisu_trim
P.phi0   = 0;  % initial roll angle
P.theta0 = 0;  % initial pitch angle
P.psi0   = 0;  % initial yaw angle
P.p0     = 0;  % initial body frame roll rate
P.q0     = 0;  % initial body frame pitch rate
P.r0     = 0;  % initial body frame yaw rate
P.delta_e0 =0;
P.delta_a0 =0;
P.delta_r0 =0;
P.delta_t0 =0;
                         

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
                         
% Aerosonde UAV Geometric Data

P.maircraft = 9.1;                   % UAV weight in kg
P.Ixx = 0.8244;                      % Moment of inertia in x-axis [kg m^2]
P.Iyy = 1.135;                       % Moment of inertia in y-axis [kg m^2]
P.Izz = 1.759;                       % Moment of inertia in z-axis [kg m^2]
P.Ixz = 0.1204;                      % Moment of inertia in xz frame [kg m^2]
P.S = 0.55;                          % Surface area of the wing in m^2
P.b = 2.8956;                        % Span of the wing in m
P.c = 0.18994;                       % Chord length of the wing in m
P.Sprop = 0.2027;                    % Disc area of the prop (in m^2)
P.R = sqrt(P.Sprop/pi);              % Prop radius
P.D = 2 * P.R;                       % Prop diameter
P.e = 0.9;                           % Oswald's efficiency factor 
P.OmegaMax_RPM = 12500;              % Max RPM 
P.Fuel_Cap_L = 4;                    % Fuel Capacity in Liters
P.Cprop = 1;                         % Aerodynamic co-efficient for propellor - Obtained from appendix E2.
P.AR = P.b^2/P.S;
P.k = 1/(pi*P.e*P.AR);

% -------------------------------------------------------------------------

% Aerodynamic Longitudinal coefficients

P.C_Lo = 0.28;
P.C_Do = 0.03;
P.C_mo = -0.02338;
P.C_Lalpha = 3.45;
P.C_Dalpha = 0.3;
P.C_malpha = -0.38;
P.C_Lq = 0;
P.C_Dq = 0;
P.C_mq = -3.6;
P.C_Ldelta_e = -0.36;
P.C_Ddelta_e = 0;
P.C_mdelta_e = -0.5;
      
% Aerodynamic Lateral Coefficients

P.C_Yo = 0;
P.C_lo = 0;
P.C_no = 0;
P.C_Ybeta = -0.98;
P.C_lbeta = -0.12;
P.C_nbeta = 0.25;
P.C_Yp = 0;
P.C_lp = -0.26;
P.C_np = 0.022;
P.C_Yr = 0;
P.C_lr = 0.14;
P.C_nr = -0.35;
P.C_Ydelta_a = 0;
P.C_ldelta_a = 0.08;
P.C_ndelta_a = 0.06;
P.C_Ydelta_r = -0.17;
P.C_ldelta_r = 0.105;
P.C_ndelta_r = -0.032;

%turn radius = 250; % turn radius in meters
%P.bank_max should be 15% of TAS(70km/h)
P.bank_max = 0.1832596; %10.5 degrees