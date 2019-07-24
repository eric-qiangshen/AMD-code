% 2016-09

% modes is a list of structs with elts:
%   name, sys, bounds
% where name is a string sys is an instance of the StateSpace class
% and bounds is a struct wit elts:
%   P_y, p_y, Q_w, q_w
% bounds is a struct with elts:
%   P_0, p_0, P_x, p_x, Q_u, q_u

function [ustar,sol] = optimal_input_noise(modes, bounds, T, epsi,NNorm)

% Infer some dimensions and stuff
N = size(modes,2);
p = size(modes(1).sys.mode(1).C,1);
m = size(modes(1).sys.mode(1).B,2);
m_u = size(bounds.Q_u,2);
m_w = size(bounds.Q_w,2);
m_v = size(bounds.Q_v,2);
m_d = m-m_w-m_u-m_v;
n = size(modes(1).sys.mode(1).A,2);
n_x = size(bounds.P_x,2);
n_y = n-n_x;
% Define optimization variables
u = sdpvar(m_u,T,'full');
s = sdpvar(p,T,N*(N-1)/2,2,'full');
a = binvar(p,T,N*(N-1)/2,2,'full');
% Building the needed matrices
Qbar_u = kron(eye(T), bounds.Q_u);
qbar_u = kron(ones(T,1), bounds.q_u);
Pbar_x = kron(eye(T*N), bounds.P_x);
pbar_x = kron(ones(T*N,1), bounds.p_x);
Pbar_y = [];
pbar_y = [];
Qbar_d=[];
qbar_d = [];
Qbar_w = [];
qbar_w = [];
Qbar_v=[];
qbar_v = [];

ftilde_x = [];
ftilde_y = [];
ftilde = [];

M_x = [];
M_y = [];
Gamma_xu = [];
Gamma_xd = [];
Gamma_xw = [];
Gamma_yu = [];
Gamma_yd = [];
Gamma_yw = [];
Gamma_u_conc = [];
Gamma_d_conc = [];
Gamma_w_conc = [];
Abar_conc = [];
E = [];
F_u=[];
F_d=[];
F_v=[];
% Loop through all modes
for i=1:N
  % Concatenate the polytope matrices
  Pbar_y = blkdiag(Pbar_y, kron(eye(T), modes(i).bounds.P_y));
  pbar_y = [pbar_y; kron(ones(T,1), modes(i).bounds.p_y)];
  Qbar_d = blkdiag(Qbar_d, kron(eye(T), modes(i).bounds.Q_d));
  qbar_d = [qbar_d; kron(ones(T,1), modes(i).bounds.q_d)];
  Qbar_w = blkdiag(Qbar_w, kron(eye(T), modes(i).bounds.Q_w));
  qbar_w = [qbar_w; kron(ones(T,1), modes(i).bounds.q_w)];
  Qbar_v = blkdiag(Qbar_v, kron(eye(T), modes(i).bounds.Q_v));
  qbar_v = [qbar_v; kron(ones(T,1), modes(i).bounds.q_v)];
  
  % Construct the matrices that multiply the initial condition to give [x; y]
  Abar_m1(:,:,i) = ones(n*(T-1), n);
  Abar(:,:,i) = ones(n*T,n);
  for k=1:(T-1)
    Abar_m1((k-1)*n+1:k*n,:,i) ...
      = modes(i).sys.mode(1).A^k;
    Abar((k-1)*n+1:k*n,:,i) ...
      = modes(i).sys.mode(1).A^k;
  end
  Abar((T-1)*n+1:T*n,:,i) = modes(i).sys.mode(1).A^T;
  % Concatenate across modes to get an i-independent Abar
  Abar_conc = blkdiag(Abar_conc, Abar(:,:,i));

  % Construct the matrices that multiply u,w,f and gives [x; y]
  % m1 stands for T 'minus' 1
   for k=1:T-1
    for kk=1:k
      Gamma_u_m1((k-1)*n+1:k*n,(kk-1)*m_u+1:kk*m_u,i) = ...
          modes(i).sys.mode(1).A^(k-kk) * modes(i).sys.mode(1).B(:,1:m_u);
      Gamma_d_m1((k-1)*n+1:k*n,(kk-1)*m_d+1:kk*m_d,i) = ...
          modes(i).sys.mode(1).A^(k-kk) * modes(i).sys.mode(1).B(:,m_u+1:m_u+m_d);
      Gamma_w_m1((k-1)*n+1:k*n,(kk-1)*m_w+1:kk*m_w,i) = ...
          modes(i).sys.mode(1).A^(k-kk) * modes(i).sys.mode(1).B(:,m_u+m_d+1:m_u+m_d+m_w);
      Gamma_u((k-1)*n+1:k*n,(kk-1)*m_u+1:kk*m_u,i) = ...
          modes(i).sys.mode(1).A^(k-kk) * modes(i).sys.mode(1).B(:,1:m_u);
      Gamma_d((k-1)*n+1:k*n,(kk-1)*m_d+1:kk*m_d,i) = ...
          modes(i).sys.mode(1).A^(k-kk) * modes(i).sys.mode(1).B(:,m_u+1:m_u+m_d);
       Gamma_w((k-1)*n+1:k*n,(kk-1)*m_w+1:kk*m_w,i) = ...
          modes(i).sys.mode(1).A^(k-kk) * modes(i).sys.mode(1).B(:,m_u+m_d+1:m_u+m_d+m_w);
      Theta_m1((k-1)*n+1:k*n,(kk-1)*n+1:kk*n,i) = modes(i).sys.mode(1).A^(k-kk);
      Theta((k-1)*n+1:k*n,(kk-1)*n+1:kk*n,i) = modes(i).sys.mode(1).A^(k-kk);
    end
  end
  
    % Then concatenate in the remaining elts for the larger matrices
  for kk=1:T
    Gamma_u((T-1)*n+1:T*n,(kk-1)*m_u+1:kk*m_u,i) = ...
        modes(i).sys.mode(1).A^(T-kk) * modes(i).sys.mode(1).B(:,1:m_u);
    Gamma_d((T-1)*n+1:T*n,(kk-1)*m_d+1:kk*m_d,i) = ...
        modes(i).sys.mode(1).A^(T-kk) * modes(i).sys.mode(1).B(:,m_u+1:m_u+m_d);
    Gamma_w((T-1)*n+1:T*n,(kk-1)*m_w+1:kk*m_w,i) = ...
        modes(i).sys.mode(1).A^(T-kk) * modes(i).sys.mode(1).B(:,1+m_u+m_d:m_u+m_d+m_w);
    Theta((T-1)*n+1:T*n,(kk-1)*n+1:kk*n,i) = modes(i).sys.mode(1).A^(T-kk);
  end


  % Build non-i-dependent, concatenated matrices
  Gamma_u_conc = [Gamma_u_conc; Gamma_u(:,:,i)];
  Gamma_d_conc = blkdiag(Gamma_d_conc, Gamma_d(:,:,i));
  Gamma_w_conc = blkdiag(Gamma_w_conc, Gamma_w(:,:,i));

  % Construct all the blockdiagonal matrices needed below
  % A_x
  A_xd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).A(1:n_x,:));
  %B_x
  B_xd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(1:n_x,:));
  B_xud(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(1:n_x,1:m_u));
  B_xdd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(1:n_x,m_u+1:m_u+m_d));
  B_xwd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(1:n_x,m_u+m_d+1:m_w+m_u+m_d));
  % A_y
  A_yd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).A(n_x+1:n,:));
  % B_y
  B_yd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(n_x+1:n,:));
  B_yud(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(n_x+1:n,1:m_u));
  B_ydd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(n_x+1:n,m_u+1:m_d+m_u));
  B_ywd(:,:,i) = kron(eye(T), modes(i).sys.mode(1).B(n_x+1:n,m_u+m_d+1:m_w+m_u+m_d));
  
  % The matrices multiplying the initial conditions giving x and y
  M_x = blkdiag(M_x,A_xd(:,:,i)*[eye(n); Abar_m1(:,:,i)]);                                             %%% modified 
  M_y = blkdiag(M_y,A_yd(:,:,i)*[eye(n); Abar_m1(:,:,i)]);                                             %%% modified 
  % The matrices multiplying u,w to give x and y
  Gamma_xu = [Gamma_xu; ...
    A_xd(:,:,i)*[zeros(n,m_u*(T-1)) zeros(n,m_u); Gamma_u_m1(:,:,i) zeros(n*(T-1),m_u)] + B_xud(:,:,i)];
  Gamma_xd = blkdiag(Gamma_xd, ...
    A_xd(:,:,i)*[zeros(n,m_d*(T-1)) zeros(n,m_d); Gamma_d_m1(:,:,i) zeros(n*(T-1),m_d)] + B_xdd(:,:,i));
  Gamma_xw = blkdiag(Gamma_xw, ...
    A_xd(:,:,i)*[zeros(n,m_w*(T-1)) zeros(n,m_w); Gamma_w_m1(:,:,i) zeros(n*(T-1),m_w)] + B_xwd(:,:,i));
  Gamma_yu = [Gamma_yu; ...
    A_yd(:,:,i)*[zeros(n,m_u*(T-1)) zeros(n,m_u); Gamma_u_m1(:,:,i) zeros(n*(T-1),m_u)] + B_yud(:,:,i)];
  Gamma_yd = blkdiag(Gamma_yd, ...
    A_yd(:,:,i)*[zeros(n,m_d*(T-1)) zeros(n,m_d); Gamma_d_m1(:,:,i) zeros(n*(T-1),m_d)] + B_ydd(:,:,i));
  Gamma_yw = blkdiag(Gamma_yw, ...
    A_yd(:,:,i)*[zeros(n,m_w*(T-1)) zeros(n,m_w); Gamma_w_m1(:,:,i) zeros(n*(T-1),m_w)] + B_ywd(:,:,i));
  % Construct the affine terms for the time-concatenated x and y states
  fbar_m1(:,i) = kron(ones(T-1,1), modes(i).sys.mode(1).f);
  fbar_x(:,i) = kron(ones(T,1), modes(i).sys.mode(1).f(1:n_x));
  fbar_y(:,i) = kron(ones(T,1), modes(i).sys.mode(1).f(n_x+1:n));
  
  ftilde_x = [ftilde_x; ...
    A_xd(:,:,i)*[zeros(n, n*(T-1)); Theta_m1(:,:,i)]*fbar_m1(:,i)+fbar_x(:,i)];
  ftilde_y = [ftilde_y; ...
    A_yd(:,:,i)*[zeros(n, n*(T-1)); Theta_m1(:,:,i)]*fbar_m1(:,i)+fbar_y(:,i)];

  % The affine term for time-concatenated [x; y]
  ftilde = [ftilde; Theta(:,:,i)*kron(ones(T,1),modes(i).sys.mode(1).f)];

  % Construct the separability-of-outputs difference-matrix, E
  for ii=(i+1):N
      E = [E; [zeros(p*T,(i-1)*n*T) kron(eye(T),modes(i).sys.mode(1).C) ...
        zeros(p*T,(ii-i-1)*n*T) -kron(eye(T),modes(ii).sys.mode(1).C) zeros(p*T,(N-ii)*n*T)]];
  end
  for ii=(i+1):N
      F_u  = [F_u; kron(eye(T), modes(i).sys.mode(1).D(:,1:m_u))- kron(eye(T), modes(ii).sys.mode(1).D(:,1:m_u))];
  end
  for ii=(i+1):N
      F_d = [F_d; [zeros(p*T,(i-1)*m_d*T) kron(eye(T),modes(i).sys.mode(1).D(:,m_u+1:m_u+m_d)) ...
        zeros(p*T,(ii-i-1)*m_d*T) -kron(eye(T),modes(ii).sys.mode(1).D(:,m_u+1:m_u+m_d)) zeros(p*T,(N-ii)*m_d*T)]];
  end
  for ii=(i+1):N
      F_v = [F_v; [zeros(p*T,(i-1)*m_v*T) kron(eye(T),modes(i).sys.mode(1).D(:,m_u+m_d+m_w+1:m_u+m_v+m_w+m_d)) ...
        zeros(p*T,(ii-i-1)*m_v*T) -kron(eye(T),modes(ii).sys.mode(1).D(:,m_u+m_d+m_w+1:m_u+m_v+m_w+m_d)) zeros(p*T,(N-ii)*m_v*T)]];
  end  
  
  
%   for ii=(i+1):N
%     E = [E; kron(eye(T), [zeros(p,(i-1)*n) modes(i).sys.mode(1).C ...
%       zeros(p,(ii-i-1)*n) -modes(ii).sys.mode(1).C zeros(p,(N-ii)*n)])];
%   end
end

% Extend E to handle both signs of output differences
Ebar = [E; -E];
F_ubar = [F_u; -F_u];
F_dbar = [F_d; -F_d];
F_vbar = [F_v; -F_v];

% Concatenate matrices further to obtain a concise optimization problem
H_x=Pbar_x*[M_x Gamma_xd Gamma_xw zeros(size(Gamma_xw,1),size(Gamma_xw,2)/m_w*m_v)];
H_y=Pbar_y*[M_y Gamma_yd Gamma_yw zeros(size(Gamma_yw,1),size(Gamma_yw,2)/m_w*m_v)]; 

h_x=pbar_x - Pbar_x*ftilde_x-Pbar_x*Gamma_xu*u(:);
h_y=pbar_y - Pbar_y*ftilde_y-Pbar_y*Gamma_yu*u(:);

H_xbar = blkdiag(kron(eye(N),bounds.P_0),Qbar_d,Qbar_w,Qbar_v);
h_xbar = [kron(ones(N,1),bounds.p_0);qbar_d;qbar_w;qbar_v]; 

G = Ebar*[Abar_conc Gamma_d_conc Gamma_w_conc zeros(size(Gamma_d_conc,1),size(Gamma_w_conc,2)/m_w*m_v)]...
    +[zeros(size(F_vbar,1),size(Abar_conc,2)),F_dbar,zeros(size(F_vbar,1),size(F_vbar,2)/m_v*m_w),F_vbar];   
g = epsi - Ebar*ftilde - (Ebar*Gamma_u_conc+F_ubar)*u(:) - s(:);

% Final matrices needed are: Qbar_u, qbar_u, Phi, phi, R, r
Phi = [H_y; H_xbar];
phi = [h_y; h_xbar];

% Final concatenations
R = [-G; H_x];
r = [-g; h_x];
%% To check whether the set of uncertaities is not empty
% xbar = sdpvar(n+(m_w*T)*N,1);
% Const = Phi*xbar<=phi;
options = sdpsettings('verbose', 1, 'solver', 'gurobi');
% sol1 = optimize(Const, [], options);
% 
% if strcmp(sol1.info, 'Infeasible problem (CPLEX-IBM)')
%     error('The feasibility set is empty')
% end

% v_max = 50/3.6; max_distance = -100;
% Dual variables for robustification
Pi = sdpvar(size(Phi,1), size(r,1),'full');
% xbar=[max_distance; v_max; max_distance; v_max;zeros((m_w*T)*N,1)];%

% Formulating the MILP

% Introduce SOS-1 constraints
sos_constraint = [];
for ij=1:(N*(N-1)/2)
  for k=1:T
    for l=1:p
      for alpha=1:2
        sos_constraint = [sos_constraint;sos1([s(l,k,ij,alpha),a(l,k,ij,alpha)],[2,1])];
      end
    end
  end
end

% Construct the constraints
% constraints = [Qbar_u*u(:)<=qbar_u;Pi>=0;transpose(Pi)*phi<=r; transpose(Phi)*Pi==transpose(R);sos_constraint];
% constraints = [Qbar_u*u(:)<=qbar_u;sos_constraint];
% constraints= [constraints;R*xbar<=r];
% constraints= [constraints;Phi*xbar<=phi];
constraints = [Qbar_u*u(:)<=qbar_u;sos_constraint];
for ij=1:(N*(N-1)/2)
  constraints = [constraints;sum(sum(sum(a(:,:,ij,:))))>=1];
end

constraints1=[constraints;Pi(:)>=0;Pi'*phi<=r;Phi'*Pi==R'];
% Objective function to minimize
%objective = u(:)'*u(:); % TODO: different obj.
if NNorm==1 || NNorm==inf
    objective = norm(u(:),NNorm);
elseif NNorm==2
    objective=u(:)'*u(:); 
elseif NNorm == 3
    objective = norm(u(:),1) + 2*norm(u(:),inf);
end
%objective = 0;
% options = sdpsettings('verbose', 1, 'solver', 'cplex');
% Solve the problem
t = tic;
sol = optimize(constraints1', objective, options);
toc(t);
% xbar2 = sdpvar(n+(m_w*T)*N,1);
% EPS = sdpvar(1,1);
% constraints2=[norm(R*xbar2-value(Pi'*phi),inf)<=EPS];
% 
% sol2 = optimize(constraints2', EPS, options);
% value(xbar2)

% Output if everything went well
if sol.problem == 0
  ustar = value(u);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
 ustar=inf;
end

end

% EOF
