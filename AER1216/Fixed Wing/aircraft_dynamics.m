function aircraft_dynamics(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

% AER1216 Fall 2023 
% Fixed Wing Project Code
%
% aircraft_dynamics.m
%
% Fixed wing simulation model file, based on the Aerosonde UAV, with code
% structure adapted from Small Unmanned Aircraft: Theory and Practice by 
% R.W. Beard and T. W. McLain. 
% 
% Inputs: 
% delta_e           elevator deflection [deg]
% delta_a           aileron deflection [deg]
% delta_r           rudder deflection [deg]
% delta_t           normalized thrust []
%
% Outputs:
% pn                inertial frame x (north) position [m]
% pe                inertial frame y (east) position [m]
% pd                inertial frame z (down) position [m]
% u                 body frame x velocity [m/s]
% v                 body frame y velocity [m/s]
% w                 body frame z velocity [m/s]
% phi               roll angle [rad]
% theta             pitch angle [rad]
% psi               yaw angle [rad]
% p                 roll rate [rad/s]
% q                 pitch rate [rad/s]
% r                 yaw rate [rad/s]
%
% Last updated: Pravin Wedage 2023-10-26

%% TA NOTE
% The main code segements you must modify are located in the derivatives
% function in this .m file. In addition, for Q7, you may need to modify the
% setup function in order to input wind into the dynamics. 
% 
% Modify other sections at your own risk. 


%
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
%
setup(block);

end 


%% Function: setup ===================================================
% Abstract:
%   Set up the basic characteristics of the S-function block such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
%
%   Required         : Yes
%   C-Mex counterpart: mdlInitializeSizes
%
function setup(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:block.NumInputPorts
    block.InputPort(i).Dimensions        = 4;
    block.InputPort(i).DatatypeID  = 0;  % double
    block.InputPort(i).Complexity  = 'Real';
    block.InputPort(i).DirectFeedthrough = false; % important to be false 
end

% Override output port properties
for i = 1:block.NumOutputPorts
    block.OutputPort(i).Dimensions       = 12;
    block.OutputPort(i).DatatypeID  = 0; % double
    block.OutputPort(i).Complexity  = 'Real';
%     block.OutputPort(i).SamplingMode = 'Sample';
end

% Register parameters
block.NumDialogPrms     = 1;
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Register multiple instances allowable
% block.SupportMultipleExecInstances = true;

% Register number of continuous states
block.NumContStates = 12;

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
% -----------------------------------------------------------------

% block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup); % discrete states only
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
block.RegBlockMethod('InitializeConditions',    @InitializeConditions);
% block.RegBlockMethod('Start',                   @Start); % Initialize Conditions is used
block.RegBlockMethod('Outputs',                 @Outputs); % Required
% block.RegBlockMethod('Update',                  @Update); % only required for discrete states
block.RegBlockMethod('Derivatives',             @Derivatives); % Required for continuous states
block.RegBlockMethod('Terminate',               @Terminate); % Required

end 


%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C-Mex counterpart: mdlSetWorkWidths
%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;

end


%% InitializeConditions:
%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C-MEX counterpart: mdlInitializeConditions
%
function InitializeConditions(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Initialize continuous states
block.ContStates.Data(1) = P.pn0; 
block.ContStates.Data(2) = P.pe0;
block.ContStates.Data(3) = P.pd0;
block.ContStates.Data(4) = P.u0;
block.ContStates.Data(5) = P.v0;
block.ContStates.Data(6) = P.w0;
block.ContStates.Data(7) = P.phi0;
block.ContStates.Data(8) = P.theta0;
block.ContStates.Data(9) = P.psi0;
block.ContStates.Data(10) = P.p0;
block.ContStates.Data(11) = P.q0;
block.ContStates.Data(12) = P.r0;

end 

%% Start:
%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this 
%                      is the place to do it.
%   Required         : No
%   C-MEX counterpart: mdlStart
%
function Start(block)

block.Dwork(1).Data = 0;

end 

%% Input Port Sampling Method:
function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = 'Sample';
  for i = 1:block.NumOutputPorts
    block.OutputPort(i).SamplingMode  = 'Sample';   
  end
end

%% Outputs:
%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C-MEX counterpart: mdlOutputs
%
function Outputs(block)

temp_mat = zeros(block.NumContStates,1); % thirteen states
for i = 1:block.NumContStates
     temp_mat(i) = block.ContStates.Data(i);
end

block.OutputPort(1).Data = temp_mat; % states

% for i = 1:block.NumOutputPorts
%     block.OutputPort(1).Data(i) = block.ContStates.Data(i);
% end

end 


%% Update:
%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C-MEX counterpart: mdlUpdate
%
function Update(block)

block.Dwork(1).Data = block.InputPort(1).Data;

end 


%% Derivatives:
%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C-MEX counterpart: mdlDerivatives
%
function Derivatives(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function


% map states and inputs
pn    = block.ContStates.Data(1);
pe    = block.ContStates.Data(2);
pd    = block.ContStates.Data(3);
u     = block.ContStates.Data(4);
v     = block.ContStates.Data(5);
w     = block.ContStates.Data(6);
phi   = block.ContStates.Data(7);
theta = block.ContStates.Data(8);
psi   = block.ContStates.Data(9);
p     = block.ContStates.Data(10);
q     = block.ContStates.Data(11);
r     = block.ContStates.Data(12);
delta_e = block.InputPort(1).Data(1)*pi/180 ; % converted inputs to radians
delta_a = block.InputPort(1).Data(2)*pi/180 ; % converted inputs to radians
delta_r = block.InputPort(1).Data(3)*pi/180 ; % converted inputs to radians
delta_t = block.InputPort(1).Data(4);


Va0=P.Va0;
g=P.g;
rho = 1.225;

% Air Data 
Va = sqrt(u^2 + v^2 + w^2);
alpha0 = atan2(w, u);
beta0 = asin(v/Va);
[~, ~, ~, rho] = atmosisa(-pd);

% Parameters of the Aerosonde UAV

ma = P.maircraft;
Ixx = P.Ixx;
Iyy = P.Iyy;
Izz = P.Izz;
Ixz = P.Ixz;
S=P.S;
b=P.b;
c=P.c;
Sprop=P.Sprop;
e = P.e;
Omega_max=P.OmegaMax_RPM; % rad/s
Fuel_Cap_L = P.Fuel_Cap_L; % L
Cprop = P.Cprop;
k = P.k;

Weight = ((ma + Fuel_Cap_L)*g);

dyn_pressure = (1/2)*rho*(Va)^2;

% Lateral Aerodynamic Coefficients 

CY0=P.C_Yo;
Cl0=P.C_lo;
Cn0=P.C_no;
CYbeta=P.C_Ybeta;
Clbeta=P.C_lbeta;
Cnbeta=P.C_nbeta;
CYp=P.C_Yp;
Clp=P.C_lp;
Cnp=P.C_np;
CYr=P.C_Yr;
Clr=P.C_lr;
Cnr=P.C_nr;
CYdelta_a=P.C_Ydelta_a;
Cldelta_a=P.C_ldelta_a;
Cndelta_a=P.C_ndelta_a;
CYdelta_r=P.C_Ydelta_r;
Cldelta_r=P.C_ldelta_r;
Cndelta_r=P.C_ndelta_r;

% Longitudinal Aerodynamic Coefficients

CL0=P.C_Lo;
CD0=P.C_Do;
Cm0=P.C_mo;
Clalpha=P.C_Lalpha;
Cdalpha=P.C_Dalpha;
Cmalpha=P.C_malpha;
Clq=P.C_Lq;
CDq=P.C_Dq;
Cmq=P.C_mq;
Cldelta_e=P.C_Ldelta_e;
CDdelta_e=P.C_Ddelta_e;
Cmdelta_e=P.C_mdelta_e;

% Computed CL and CD values.
CL = CL0 + Clalpha * alpha0;
CD = CD0 + Cdalpha * alpha0;

% Compute inertial constants (Capital Gamma)

K = ((Ixx*Izz)-(Ixz^2));
k1 = (Ixz*(Ixx-Iyy+Izz))/K;
k2 = ((Izz*(Izz-Iyy))+(Ixz^2))/K;
k3 = Izz/K;
k4 = Ixz/K;
k5 = (Izz-Ixx)/Iyy;
k6 = Ixz/Iyy;
k7 = (((Ixx-Iyy)*Ixx)+(Ixz^2))/K;
k8 = Ixx/K;

% Non-Linear Dynamics

% Rotation matrix
C_BE = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
        sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);
        cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)];

%SB Inverse Matrix
S_B_inv = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
           0, cos(phi), -sin(phi);
           0, sin(phi)*sec(theta), cos(phi)*sec(theta)];


% Aerodynamic forces and moments
F_lift = dyn_pressure*S*(CL + (Clq*(c/(2*Va))*q) + (Cldelta_e*delta_e)); % lift force at alpha=0
F_drag = dyn_pressure*S*(CD + (CDq*(c/(2*Va))*q) + (CDdelta_e*delta_e)); % drag force
Fy = dyn_pressure*S*(CY0 + (CYbeta*beta0) + (CYp*(b/(2*Va))*p) + (CYr*(b/(2*Va))*r) + (CYdelta_a*delta_a) + (CYdelta_r*delta_r)); % sideforce

m = (dyn_pressure*S*c*(Cm0 + (Cmalpha*alpha0) + (Cmq*(c/(2*Va))*q) + (Cmdelta_e*delta_e))); % pitching moment
l = (dyn_pressure*S*b*(Cl0 + (Clbeta*beta0) + (Clp*(b/(2*Va))*p) + (Clr*(b/(2*Va))*r) + (Cldelta_a*delta_a) + (Cldelta_r*delta_r))); % rolling moment 
n = (dyn_pressure*S*b*(Cn0 + (Cnbeta*beta0) + (Cnp*(b/(2*Va))*p) + (Cnr*(b/(2*Va))*r) + (Cndelta_a*delta_a) + (Cndelta_r*delta_r))); % yawing moment

% Compute the gravitational forces here
Fg_x = -Weight*sin(theta);
Fg_y = Weight*cos(theta)*sin(phi);
Fg_z = Weight*cos(theta)*cos(phi);

% Propulsion Force
%Fprop = (0.5*rho*Sprop*Cprop*(((80*delta_t)^2) - (Va)^2));
%Fprop = Ct*1.225*(n^2)*(0.5080^4);
% Fprop = T;
Fprop = 177.04*(delta_t);

% Total forces and moments (body frame)
Stability2BodyFrame = [ cos(alpha0), -sin(alpha0);
                        sin(alpha0), cos(alpha0)];

F_long = Stability2BodyFrame * [-F_drag; -F_lift];

X = Fg_x + F_long(1) + Fprop;
Y = Fg_y + Fy;
Z = Fg_z + F_long(2);
L = l;
M = m;
N = n;

% State derivatives
% The full aircraft dynamics model is computed here
pdot = C_BE' * [u; v; w];
pndot = pdot(1);
pedot = pdot(2);
pddot = pdot(3);

udot = (X/(ma+Fuel_Cap_L)) + (r*v) - (q*w);
vdot = (Y/(ma+Fuel_Cap_L)) + (p*w) - (r*u);
wdot = (Z/(ma+Fuel_Cap_L)) + (q*u) - (p*v);

euler_angles_dot = S_B_inv * [p; q; r];
phidot = euler_angles_dot(1);
thetadot = euler_angles_dot(2);
psidot = euler_angles_dot(3);

angular_rates = [((k1*p*q) - (k2*q*r) + (k3*l) + (k4*n));
                 ((k5*p*r) - (k6*(p^2-r^2)) + ((1/Iyy)*m));
                 ((k7*p*q) - (k1*q*r) + (k4*l) + (k8*n))];

pdot = angular_rates(1);
qdot = angular_rates(2);
rdot = angular_rates(3);


% map derivatives
block.Derivatives.Data(1) = pndot;
block.Derivatives.Data(2) = pedot;
block.Derivatives.Data(3) = pddot;
block.Derivatives.Data(4) = udot;
block.Derivatives.Data(5) = vdot;
block.Derivatives.Data(6) = wdot;
block.Derivatives.Data(7) = phidot;
block.Derivatives.Data(8) = thetadot;
block.Derivatives.Data(9) = psidot;
block.Derivatives.Data(10)= pdot;
block.Derivatives.Data(11)= qdot;
block.Derivatives.Data(12)= rdot;

end 


%% Terminate:
%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C-MEX counterpart: mdlTerminate
%
function Terminate(block)

end 

