% Here is the code of active model discrimination for intersection crossing
% scenario
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
T_hor = 5; % length of horizon
epsi = 0.25; % separability threshold
NNorm = 3; % specify the norm in objective function

% Bounds
v_max_x = 9; % 30*1.6/3.6; %30mile/h
v_max_y = 9; % 30*1.6/3.6; %30mile/h
v_min_x = 6; % 17*1.6/3.6; % 20mile/h
v_min_y = 6; % 17*1.6/3.6; % 20mile/h
u_max = 100/7/3.6; % 0-100km/h in 7s
u_min = -0.8*9.81; % Totota Corolla's braking power 0.8*g

% Initial condition bounds
P_0 = [eye(4); -eye(4)]; 
p_0 = [-15;  v_max_x; -15;  v_max_y; 
        18; -v_min_x;  18; -v_min_y];

% Ego car's state bounds
P_x = [0 1; 0 -1]; 
p_x = [v_max_x; 0];
% Ego car's input bounds
Q_u = [1; -1];
q_u = [u_max; -u_min];
% The process noise
w_max = 1;
w_min = -1;
Q_w = [1 0; -1 0; 0 1; 0 -1];
q_w = 0.01*[w_max; -w_min; w_max; -w_min];
% The measurement noise
vv_max = 1;
vv_min = -1;
Q_v = [1; -1];
q_v = 0.01*[vv_max; -vv_min];
bounds = struct('P_0', P_0, 'p_0', p_0, 'P_x', P_x, 'p_x', p_x, 'Q_u', Q_u, ...
                'q_u', q_u, 'Q_w', Q_w, 'q_w', q_w, 'Q_v', Q_v, 'q_v', q_v);
% The othercar's input bound
Q_d = Q_u;
q_d = q_u;

%% Intentions (modes/models)
% Inattentive driver
% u_y = dI, where dI is input uncertaity satisfying a polyhedral set
AI(:,:,1) = [1 T 0 0; ...
             0 1 0 0; ...
             0 0 1 T; ...
             0 0 0 1];
BI(:,:,1) = [0 0 0 0 0;...
             T 0 T 0 0;...
             0 0 0 0 0;...
             0 T 0 T 0];
CI(:,:,1) = [0 0 0 1];
DI(:,:,1) = [0 0 0 0 1];
fI(:,:,1) = zeros(4,1);
sysI = StateSpace(AI, BI, CI, DI, fI, []);
% other car's state constaint
% P_yI = zeros(2,2);
% p_yI = zeros(2,1);
P_yI = [0 1; 0 -1];
p_yI = [9; -6];
Q_dI = Q_d;
q_dI = 0.1*q_d;
Q_wI = Q_w;
q_wI = q_w;
Q_vI = Q_v;
q_vI = q_v;
% The structs have the exact same form for all the models (important!),
% can easily be changed to classes later to enforce that.
bounds_I = struct('P_y', P_yI, 'p_y', p_yI, 'Q_d', Q_dI, 'q_d', q_dI, 'Q_w', Q_wI, ...
                  'q_w', q_wI, 'Q_v', Q_vI, 'q_v', q_vI);
struct_I = struct('name', 'I', 'sys', sysI, 'bounds', bounds_I);

% Malicious intention
% u_y = K_pM*(x-y) + K_dM*(v_x-x_y) + dM, 
K_dM = 3.5;
K_pM = 1;
AM(:,:,1) = [1      T       0       0; ...
             0      1       0       0; ...
             0      0       1       T; ...
             T*K_pM T*K_dM  -T*K_pM 1-T*K_dM];
BM = BI;
CM = CI;
DM = DI;
fM = fI;
sysM = StateSpace(AM, BM, CM, DM, fM, []);
P_yM = zeros(2,2); 
p_yM = zeros(2,1); 
% P_yM = [0 1; 0 -1];
% p_yM = [9; -6];
Q_dM = Q_d;
q_dM = 0.05*q_u;
Q_wM = Q_w;
q_wM = q_w;
Q_vM = Q_v;
q_vM = q_v;
bounds_M = struct('P_y', P_yM, 'p_y', p_yM, 'Q_d', Q_dM, 'q_d', q_dM, ...
                  'Q_w', Q_wM, 'q_w', q_wM, 'Q_v', Q_vM, 'q_v', q_vM);
struct_M = struct('name', 'M', 'sys', sysM, 'bounds', bounds_M);

% Cautious intention
% u_y = -K_pC*y - K_dC*v_y + dC 
K_dC = 4.75;
K_pC = 1.5;
AC(:,:,1) = [1 T 0          0; ...
             0 1 0          0; ...
             0 0 1          T; ...
             0 0 -T*K_pC    1-T*K_dC];
BC=BI;
CC=CI;
DC=DI;
fC=fI;
sysC = StateSpace(AC, BC, CC, DC, fC, []);
P_yC = zeros(2,2); 
p_yC = zeros(2,1); 
% P_yC = [0 1; 0 -1];
% p_yC = [9; -3]; 
Q_dC = Q_d;
q_dC = 0.05*q_u;
Q_wC = Q_w;
q_wC = q_w;
Q_vC = Q_v;
q_vC = q_v;
bounds_C = struct('P_y', P_yC, 'p_y', p_yC, 'Q_d', Q_dC, 'q_d', q_dC, ...
                  'Q_w', Q_wC, 'q_w', q_wC, 'Q_v', Q_vC, 'q_v', q_vC);
struct_C = struct('name', 'C', 'sys', sysC, 'bounds', bounds_C);

modes = [struct_M struct_C struct_I];

t1 = tic;

%% Exact method:
%  OUTER problem with ego repsonsibility: Robust Optimization; 
%  INNER problem: KKT

[u,sol] = Exact_get_u_respon(modes, bounds, T_hor, epsi, NNorm);

%% Conservative method: Robust Otpimization
%  Reference: Harirchi, F., Yong, S. Z., Jacobsen E., and Ozay, N. (2017)
%  Active model discrimination with applications to fraud detection in
%  smart buildings, In IFAC World Congress, Toulouse, France.

% [u,sol] = optimal_input_noise(modes, bounds, T_hor, epsi, NNorm);

%% Save results
if (NNorm == 1)
    u_int_1_time=toc(t1)
    u_int_1=u
	save('u_int_1_resp','u_int_1','u_int_1_time')
elseif (NNorm == 2)
    u_int_2_time=toc(t1)
    u_int_2=u
    save('u_int_2_resp','u_int_2','u_int_2_time')
elseif (NNorm == inf)
    u_int_inf_time=toc(t1)
    u_int_inf=u
    save('u_int_inf_resp','u_int_inf','u_int_inf_time')
elseif (NNorm == 3) % 1-norm + 2*inf-norm
    u_int_1inf_time=toc(t1)
    u_int_1inf=u
    save('u_int_1inf_resp','u_int_1inf','u_int_1inf_time')
end


