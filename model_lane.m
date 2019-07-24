% Here is the code of active model discrimination for lane changing scenario  
% 
% Reference: Ding, Y., Harirchi F., Yong, S. Z., Jacobsen, E., and Ozay, N. (2018).
% Optimal input design for affine model discrimination with applications in intention-aware vehicles, 
% in Proceedings of the 9th ACM/IEEE International Conference on
% Cyber-Physical Systems, PP. 297-307, arXiv:1702.01112
%
% 07/22/2019 

clc; clear all; close all;

%% Setup
% Time step (s)
T = 0.3; % sampling time
T_hor = 4; % length of horizon
epsi = 0.25; % separability threshold
NNorm = 1; % specify the norm in objective function

% Bounds on ego car's state
vx_e_max = 35; % 
vx_e_min = 27;
y_e_max = 2;
y_e_min = 0.5;
vy_e_max = 0;
vy_e_min = -0.35;
% Bounds on other car's state
vx_o_max = 35;
vx_o_min = 27;

u_x_max = 100/7/3.6; % 0-100km/h in 7s
u_x_min = -0.8*9.81; % Totota Corolla's braking power 0.8*g
% Initial condition bounds
P_0 = [1  0   0     0   0;...
       -1 0   0     0   0;...
       0  0   0     1   0;...
       0  0   0     -1  0;...
       0  1   0     0   0;...
       0  -1  0     0   0;...
       0  0   1     0   0;...
       0  0   -1    0   0;...
       0  0   0     0   1;...
       0  0   0     0   -1];
p_0 =[0;    0;...
      -7;   12;...
      32;   -30;...
      1.8;  -1.1;...
      32;   -30];
% Ego car's state bounds
P_x =[0 1   0;...
      0 -1  0;...
      0 0   1;...
      0 0   -1];
p_x = [vx_e_max; -vx_e_min; y_e_max; -y_e_min]; 
% Ego car's input bounds
Q_u = [1 0; -1 0; 0 1; 0 -1];
q_u = [u_x_max; -u_x_min; vy_e_max; -vy_e_min];
% The process noise
w_max = 1;
w_min = -1;
Q_w = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
q_w = 0.01*[w_max; -w_min;w_max; -w_min;w_max; -w_min];
% The measurement noise
vv_max = 1;
vv_min = -1;
Q_v = [1; -1];
q_v = 0.01*[vv_max; -vv_min];
% Structure "bounds" stores initial constraint, ego car state
% constraint, input constraint, process noise constraint, output noise
% constraint
bounds = struct('P_0', P_0, 'p_0', p_0, 'P_x', P_x, 'p_x', p_x, 'Q_u', Q_u, 'q_u', q_u, ...
                'Q_w', Q_w, 'q_w', q_w, 'Q_v', Q_v, 'q_v', q_v);
% The other car's input bounds (same as the ego car)
Q_d = [1;-1];
q_d = [u_x_max;-u_x_min];
%% Intentions (modes / models)
% Inattentive model
% u_o = dI, where dI is input uncertaity satisfying a polyhedral set
AI = [1 T 0 0 0; ...
      0 1 0 0 0;...
      0 0 1 0 0;...
      0 0 0 1 T;...
      0 0 0 0 1];
BI = [0 0 0 0 0 0 0;...
      T 0 0 T 0 0 0;...
      0 T 0 0 T 0 0;...
      0 0 0 0 0 0 0;...
      0 0 T 0 0 T 0];
CI = [0 0 0 0 1]; 
DI= [0 0 0 0 0 0 1]; 
fI = zeros(5,1);
sysI = StateSpace(AI, BI, CI, DI, fI, []);
% other car's state constaint
% P_yI = [0 1;...
%         0 -1];
% p_yI = [vx_o_max; -vx_o_min];
P_yI = zeros(2,2);
p_yI = zeros(2,1);
Q_dI = Q_d;
q_dI = 0.1*q_d;
Q_wI = Q_w;
q_wI = q_w;
Q_vI = Q_v;
q_vI = q_v;
% The structs have the exact same form for all the models
% Structure "bounds_I" stores cautious car's state constraint, uncertainty
% constraint, process noise constraint, output noise constraint
bounds_I = struct('P_y', P_yI, 'p_y', p_yI, 'Q_d', Q_dI, 'q_d', q_dI, 'Q_w', Q_wI, ...
                  'q_w', q_wI, 'Q_v', Q_vI, 'q_v', q_vI);
struct_I = struct('name', 'I', 'sys', sysI, 'bounds', bounds_I);

% Malicious driver
% u_o = K_dM*(vx_e-vx_o) + L_pM*(y_bar-y_e) - L_dM*vy_e  + dM, 
% dM is input uncertaity satisfying a polyhedral set
K_dM = 1.1;
L_dM = 8.7;
L_pM = 2.0;
y_bar = 2;
AM = [1     T        0        0   0; ...
      0     1        0        0   0;...
      0     0        1        0   0;...
      0     0        0        1   T;...
      0     T*K_dM   -T*L_pM  0   1-T*K_dM];
BM = [0   0        0   0 0 0 0;...
      T   0        0   T 0 0 0;...
      0   T        0   0 T 0 0;...
      0   0        0   0 0 0 0;...
      0   -T*L_dM  T   0 0 T 0];
CM = CI;
DM = DI;
fM = [0; 0; 0; 0; L_pM*T*y_bar];
sysM = StateSpace(AM, BM, CM, DM, fM, []);
P_yM = zeros(2,2);
p_yM = zeros(2,1);
% P_yM = [0 1;...
%         0 -1];
% p_yM = [vx_o_max; -vx_o_min];
Q_dM = Q_d;
q_dM = 0.05*q_d;
Q_wM = Q_w;
q_wM = q_w;
Q_vM = Q_v;
q_vM = q_v;
% Structure "bounds_M" stores malicious car's state constraint, uncertainty
% constraint, process noise constraint, output noise constraint
bounds_M = struct('P_y', P_yM, 'p_y', p_yM, 'Q_d', Q_dM, 'q_d', q_dM,'Q_w', Q_wM, 'q_w', q_wM,'Q_v', Q_vM, 'q_v', q_vM);
struct_M = struct('name', 'M', 'sys', sysM, 'bounds', bounds_M);

% Cautious model
% u_o = -K_dC*(vx_e-vx_o) - L_pC*(y_bar-y_e) + L_dC*vy_e + dC
K_dC = 0.9;
L_dC = 8.9;
L_pC = 2.5;
AC =[1  T           0         0   0; ...
     0  1           0         0   0;...
     0  0           1         0   0;...
     0  0           0         1   T;...
     0  -T*K_dC     T*L_pC    0   1+T*K_dC];
BC =[0 0 0 0 0 0 0;...
     T 0        0   T 0 0 0;...
     0 T        0   0 T 0 0;...
     0 0        0   0 0 0 0;...
     0 T*L_dC   T   0 0 T 0];
CC = CI;
DC = DI;
fC = [0; 0; 0; 0; -L_pC*T*y_bar];
sysC = StateSpace(AC, BC, CC, DC, fC, []);
P_yC = zeros(2,2);
p_yC = zeros(2,1);
% P_yC = [0 1;...
%         0 -1];
% p_yC = [vx_o_max; -vx_o_min];
Q_dC = Q_d;
q_dC = 0.05*q_d;
Q_wC = Q_w;
q_wC = q_w;
Q_vC = Q_v;
q_vC = q_v;
% Structure "bounds_C" stores cautious car's state constraint, uncertainty
% constraint, process noise constraint, output noise constraint
bounds_C = struct('P_y', P_yC, 'p_y', p_yC, 'Q_d', Q_dC, 'q_d', q_dC, ...
                  'Q_w', Q_wC, 'q_w', q_wC, 'Q_v', Q_vC, 'q_v', q_vC);
struct_C = struct('name', 'C', 'sys', sysC, 'bounds', bounds_C);

% Structure "modes" stores constriants on different intentions
modes = [struct_C struct_I struct_M];

t1 = tic;

%% Exact method:
%  OUTER problem with ego repsonsibility: Robust Optimization; 
%  INNER problem: KKT

% [u,sol] = Exact_get_u_respon(modes, bounds, T_hor, epsi, NNorm);

%% Conservative method: Robust Otpimization
%  Reference: Harirchi, F., Yong, S. Z., Jacobsen E., and Ozay, N. (2017)
%  Active model discrimination with applications to fraud detection in
%  smart buildings, In IFAC World Congress, Toulouse, France.

[u,sol] = optimal_input_noise(modes, bounds, T_hor, epsi, NNorm);

%% Save results
if (NNorm == 1)
    u_lane_1_time = toc(t1)
    u_lane_1 = u
	save('u_lane_1_resp', 'u_lane_1', 'u_lane_1_time')
elseif (NNorm == 2)
    u_lane_2_time = toc(t1)
    u_lane_2 = u
    save('u_lane_2_resp', 'u_lane_2', 'u_lane_2_time')
elseif (NNorm == inf)
    u_lane_inf_time = toc(t1)
    u_lane_inf = u
    save('u_lane_inf_resp','u_lane_inf','u_lane_inf_time')
elseif (NNorm == 3)
    u_lane_1inf_time = toc(t1)
    u_lane_1inf = u
    save('u_lane_1inf_resp', 'u_lane_1inf', 'u_lane_1inf_time')
end
