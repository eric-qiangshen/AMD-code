% Computing the optimal separating input for active model discrimination
%
% Reference: Ding, Y., Harirchi F., Yong, S. Z., Jacobsen, E., and Ozay, N. (2018).
% Optimal input design for affine model discrimination with applications in intention-aware vehicles, 
% in Proceedings of the 9th ACM/IEEE International Conference on
% Cyber-Physical Systems, PP. 297-307, arXiv:1702.01112
%
% 07/22/2019 

function [ustar, sol] = Exact_get_u_respon(modes, bounds, T, epsi, NNorm)

% Infer some dimensions 
N = size(modes, 2); % No. of intentions
p = size(modes(1).sys.mode(1).C, 1); % Dim of output
m = size(modes(1).sys.mode(1).B, 2); % Dim of input
m_u = size(bounds.Q_u, 2); % Dim of controlled input
m_w = size(bounds.Q_w, 2); % Dim of process noise
m_v = size(bounds.Q_v, 2); % Dim of output noise
m_d = m-m_w-m_u-m_v; % Dim of uncontrolled input
n = size(modes(1).sys.mode(1).A, 2); % Dim of total state
n_x = size(bounds.P_x, 2); % Dim of controlled state
n_y = n-n_x; % Dim of uncontrolled state

% Define optimization variables
u = sdpvar(m_u, T, 'full');

%% Building the needed matrices
% Constraints on input and controlled state
Qbar_u = kron(eye(T), bounds.Q_u);
qbar_u = kron(ones(T,1), bounds.q_u);
Pbar_x = kron(eye(T*2), bounds.P_x);
pbar_x = kron(ones(T*2,1), bounds.p_x);

% Loop through all modes to find concatenated constraints
q = 0; % modes pair index

for i = 1:N-1
    % Construct the matrices that multiply the initial condition in concatenated dynamics 
    % m1 stands for T 'minus' 1
    Abar_m1(:,:,i) = ones(n*(T-1), n);
    Abar(:,:,i) = ones(n*T, n);
    for k = 1:(T-1)
        Abar_m1((k-1)*n+1:k*n,:,i) ...
            = modes(i).sys.mode(1).A^(k);
        Abar((k-1)*n+1:k*n,:,i) ...
            = modes(i).sys.mode(1).A^(k);
    end
    Abar((T-1)*n+1:T*n,:,i) = modes(i).sys.mode(1).A^T;
    
    % First loop to T-1 for the smaller matrices
    for k = 1:T-1
        for kk = 1:k
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
    for kk = 1:T
        Gamma_u((T-1)*n+1:T*n,(kk-1)*m_u+1:kk*m_u,i) = ...
            modes(i).sys.mode(1).A^(T-kk) * modes(i).sys.mode(1).B(:,1:m_u);
        Gamma_d((T-1)*n+1:T*n,(kk-1)*m_d+1:kk*m_d,i) = ...
            modes(i).sys.mode(1).A^(T-kk) * modes(i).sys.mode(1).B(:,m_u+1:m_u+m_d);
        Gamma_w((T-1)*n+1:T*n,(kk-1)*m_w+1:kk*m_w,i) = ...
            modes(i).sys.mode(1).A^(T-kk) * modes(i).sys.mode(1).B(:,1+m_u+m_d:m_u+m_d+m_w);
        Theta((T-1)*n+1:T*n,(kk-1)*n+1:kk*n,i) = modes(i).sys.mode(1).A^(T-kk);
    end
    
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
    
    % Construct the affine terms for the time-concatenated x and y states
    fbar_m1(:,i) = kron(ones(T-1,1), modes(i).sys.mode(1).f);
    fbar_x(:,i) = kron(ones(T,1), modes(i).sys.mode(1).f(1:n_x));
    fbar_y(:,i) = kron(ones(T,1), modes(i).sys.mode(1).f(n_x+1:n));
    
    for j = i+1:N
        q = q+1; % pair index; q=1 ==> (1,2), q=2 ==> (1,3), q=1 ==> (2,3),
        % Concatenate the polytope matrices
        Pbar_y{q} = blkdiag(kron(eye(T), modes(i).bounds.P_y), kron(eye(T), modes(j).bounds.P_y));
        pbar_y{q} = [kron(ones(T,1), modes(i).bounds.p_y); kron(ones(T,1), modes(j).bounds.p_y)];
        
        Qbar_d{q} = blkdiag(kron(eye(T), modes(i).bounds.Q_d),kron(eye(T), modes(j).bounds.Q_d));
        qbar_d{q} = [kron(ones(T,1), modes(i).bounds.q_d);kron(ones(T,1), modes(j).bounds.q_d)];
        
        Qbar_w{q} = blkdiag(kron(eye(T), modes(i).bounds.Q_w),kron(eye(T), modes(j).bounds.Q_w));
        qbar_w{q} = [kron(ones(T,1), modes(i).bounds.q_w);kron(ones(T,1), modes(j).bounds.q_w)];
        
        Qbar_v{q} = blkdiag(kron(eye(T), modes(i).bounds.Q_v),kron(eye(T), modes(j).bounds.Q_v));
        qbar_v{q} = [kron(ones(T,1), modes(i).bounds.q_v);kron(ones(T,1), modes(j).bounds.q_v)];
        
        % Construct the matrices that multiply the initial condition to give [x; y]
        Abar_m1(:,:,j) = ones(n*(T-1),n);
        Abar(:,:,j) = ones(n*T,n);
        for k = 1:(T-1)
            Abar_m1((k-1)*n+1:k*n,:,j) ...
                = modes(j).sys.mode(1).A^(k);
            Abar((k-1)*n+1:k*n,:,j) ...
                = modes(j).sys.mode(1).A^(k);
        end
        Abar((T-1)*n+1:T*n,:,j) = modes(j).sys.mode(1).A^T;
        % Concatenate across modes to get an i-independent Abar
        Abar_conc{q} = blkdiag(Abar(:,:,i),Abar(:,:,j));
        
        % Construct the matrices that multiply u,d, w,f and gives [x; y]
        % m1 stands for T 'minus' 1
        % First loop to T-1 for the smaller matrices
        for k = 1:T-1
            for kk = 1:k
                Gamma_u_m1((k-1)*n+1:k*n,(kk-1)*m_u+1:kk*m_u,j) = ...
                    modes(j).sys.mode(1).A^(k-kk) * modes(j).sys.mode(1).B(:,1:m_u);
                Gamma_d_m1((k-1)*n+1:k*n,(kk-1)*m_d+1:kk*m_d,j) = ...
                    modes(j).sys.mode(1).A^(k-kk) * modes(j).sys.mode(1).B(:,m_u+1:m_u+m_d);
                Gamma_w_m1((k-1)*n+1:k*n,(kk-1)*m_w+1:kk*m_w,j) = ...
                    modes(j).sys.mode(1).A^(k-kk) * modes(j).sys.mode(1).B(:,m_u+m_d+1:m_w+m_u+m_d);
                Gamma_u((k-1)*n+1:k*n,(kk-1)*m_u+1:kk*m_u,j) = ...
                    modes(j).sys.mode(1).A^(k-kk) * modes(j).sys.mode(1).B(:,1:m_u);
                Gamma_d((k-1)*n+1:k*n,(kk-1)*m_d+1:kk*m_d,j) = ...
                    modes(j).sys.mode(1).A^(k-kk) * modes(j).sys.mode(1).B(:,m_u+1:m_u+m_d);
                Gamma_w((k-1)*n+1:k*n,(kk-1)*m_w+1:kk*m_w,j) = ...
                    modes(j).sys.mode(1).A^(k-kk) * modes(j).sys.mode(1).B(:,m_u+m_d+1:m_w+m_u+m_d);
                Theta_m1((k-1)*n+1:k*n,(kk-1)*n+1:kk*n,j) = modes(j).sys.mode(1).A^(k-kk);
                Theta((k-1)*n+1:k*n,(kk-1)*n+1:kk*n,j) = modes(j).sys.mode(1).A^(k-kk);
            end
        end
        
        % Then concatenate in the remaining elts for the larger matrices
        for kk = 1:T
            Gamma_u((T-1)*n+1:T*n,(kk-1)*m_u+1:kk*m_u,j) = ...
                modes(j).sys.mode(1).A^(T-kk) * modes(j).sys.mode(1).B(:,1:m_u);
            Gamma_d((T-1)*n+1:T*n,(kk-1)*m_d+1:kk*m_d,j) = ...
                modes(j).sys.mode(1).A^(T-kk) * modes(j).sys.mode(1).B(:,m_u+1:m_u+m_d);
            Gamma_w((T-1)*n+1:T*n,(kk-1)*m_w+1:kk*m_w,j) = ...
                modes(j).sys.mode(1).A^(T-kk) * modes(j).sys.mode(1).B(:,m_u+m_d+1:m_w+m_u+m_d);
            Theta((T-1)*n+1:T*n,(kk-1)*n+1:kk*n,j) = modes(j).sys.mode(1).A^(T-kk);
        end
        
        % Build non-i-dependent, concatenated matrices
        Gamma_u_conc{q} = [Gamma_u(:,:,i);Gamma_u(:,:,j)];
        Gamma_d_conc{q} = blkdiag(Gamma_d(:,:,i),Gamma_d(:,:,j));
        Gamma_w_conc{q} = blkdiag(Gamma_w(:,:,i),Gamma_w(:,:,j));
        
        % Construct all the blockdiagonal matrices needed below
        A_xd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).A(1:n_x,:));
        B_xd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(1:n_x,:));
        B_xud(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(1:n_x,1:m_u));
        B_xdd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(1:n_x,m_u+1:m_u+m_d));
        B_xwd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(1:n_x,m_u+m_d+1:m_u+m_d+m_w));
        A_yd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).A(n_x+1:n,:));
        B_yd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(n_x+1:n,:));
        B_yud(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(n_x+1:n,1:m_u));
        B_ydd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(n_x+1:n,m_u+1:m_u+m_d));
        B_ywd(:,:,j) = kron(eye(T), modes(j).sys.mode(1).B(n_x+1:n,m_u+m_d+1:m_u+m_d+m_w));
        
        % The matrices multiplying the initial conditions giving x and y
        M_x{q} = blkdiag(A_xd(:,:,i)*[eye(n); Abar_m1(:,:,i)],A_xd(:,:,j)*[eye(n); Abar_m1(:,:,j)]);
        M_y{q} = blkdiag(A_yd(:,:,i)*[eye(n); Abar_m1(:,:,i)],A_yd(:,:,j)*[eye(n); Abar_m1(:,:,j)]);
        
        % The matrices multiplying u,d,w to give x and y
        Gamma_xu{q} = [A_xd(:,:,i)*[zeros(n,m_u*(T-1)) zeros(n,m_u); Gamma_u_m1(:,:,i) zeros(n*(T-1),m_u)] + B_xud(:,:,i);...
            A_xd(:,:,j)*[zeros(n,m_u*(T-1)) zeros(n,m_u); Gamma_u_m1(:,:,j) zeros(n*(T-1),m_u)] + B_xud(:,:,j)];
        
        Gamma_xd{q} = blkdiag(A_xd(:,:,i)*[zeros(n,m_d*(T-1)) zeros(n,m_d); Gamma_d_m1(:,:,i) zeros(n*(T-1),m_d)] + B_xdd(:,:,i),...
            A_xd(:,:,j)*[zeros(n,m_d*(T-1)) zeros(n,m_d); Gamma_d_m1(:,:,j) zeros(n*(T-1),m_d)] + B_xdd(:,:,j));
        
        Gamma_xw{q} = blkdiag(A_xd(:,:,i)*[zeros(n,m_w*(T-1)) zeros(n,m_w); Gamma_w_m1(:,:,i) zeros(n*(T-1),m_w)] + B_xwd(:,:,i),...
            A_xd(:,:,j)*[zeros(n,m_w*(T-1)) zeros(n,m_w); Gamma_w_m1(:,:,j) zeros(n*(T-1),m_w)] + B_xwd(:,:,j));
        
        Gamma_yu{q} = [A_yd(:,:,i)*[zeros(n,m_u*(T-1)) zeros(n,m_u); Gamma_u_m1(:,:,i) zeros(n*(T-1),m_u)] + B_yud(:,:,i);...
            A_yd(:,:,j)*[zeros(n,m_u*(T-1)) zeros(n,m_u); Gamma_u_m1(:,:,j) zeros(n*(T-1),m_u)] + B_yud(:,:,j)];
        
        Gamma_yd{q} = blkdiag(A_yd(:,:,i)*[zeros(n,m_d*(T-1)) zeros(n,m_d); Gamma_d_m1(:,:,i) zeros(n*(T-1),m_d)] + B_ydd(:,:,i),...
            A_yd(:,:,j)*[zeros(n,m_d*(T-1)) zeros(n,m_d); Gamma_d_m1(:,:,j) zeros(n*(T-1),m_d)] + B_ydd(:,:,j));
        
        Gamma_yw{q} = blkdiag(A_yd(:,:,i)*[zeros(n,m_w*(T-1)) zeros(n,m_w); Gamma_w_m1(:,:,i) zeros(n*(T-1),m_w)] + B_ywd(:,:,i),...
            A_yd(:,:,j)*[zeros(n,m_w*(T-1)) zeros(n,m_w); Gamma_w_m1(:,:,j) zeros(n*(T-1),m_w)] + B_ywd(:,:,j));
 
        
        % Construct the affine terms for the time-concatenated x and y states
        fbar_m1(:,j) = kron(ones(T-1,1), modes(j).sys.mode(1).f);
        fbar_x(:,j) = kron(ones(T,1), modes(j).sys.mode(1).f(1:n_x));
        fbar_y(:,j) = kron(ones(T,1), modes(j).sys.mode(1).f(n_x+1:n));
        
        ftilde_x{q} = [A_xd(:,:,i)*[zeros(n, n*(T-1)); Theta_m1(:,:,i)]*fbar_m1(:,i)+fbar_x(:,i);...
            A_xd(:,:,j)*[zeros(n, n*(T-1)); Theta_m1(:,:,j)]*fbar_m1(:,j)+fbar_x(:,j)];
        ftilde_y{q}  = [A_yd(:,:,i)*[zeros(n, n*(T-1)); Theta_m1(:,:,i)]*fbar_m1(:,i)+fbar_y(:,i);...
            A_yd(:,:,j)*[zeros(n, n*(T-1)); Theta_m1(:,:,j)]*fbar_m1(:,j)+fbar_y(:,j)];
        
        % The affine term for time-concatenated [x; y]
        ftilde{q}  = [Theta(:,:,i)*kron(ones(T,1),modes(i).sys.mode(1).f);
            Theta(:,:,j)*kron(ones(T,1),modes(j).sys.mode(1).f)];
        
        % Construct the separability-of-outputs difference-matrix, E
        E{q} = [kron(eye((T)), modes(i).sys.mode(1).C), -kron(eye((T)), modes(j).sys.mode(1).C); ...
                -kron(eye((T)), modes(i).sys.mode(1).C), kron(eye((T)), modes(j).sys.mode(1).C)];
        
        F_u{q} = [kron(eye((T)), modes(i).sys.mode(1).D(:,1:m_u)) - kron(eye((T)), modes(j).sys.mode(1).D(:,1:m_u)); ...
                  -kron(eye((T)), modes(i).sys.mode(1).D(:,1:m_u)) + kron(eye((T)), modes(j).sys.mode(1).D(:,1:m_u))];
        
        F_d{q} = [kron(eye((T)), modes(i).sys.mode(1).D(:,m_u+1:m_u+m_d)), -kron(eye((T)), modes(j).sys.mode(1).D(:,m_u+1:m_u+m_d)); ...
                  -kron(eye((T)), modes(i).sys.mode(1).D(:,m_u+1:m_u+m_d)), kron(eye((T)), modes(j).sys.mode(1).D(:,m_u+1:m_u+m_d))];
        
        F_v{q} = [kron(eye((T)), modes(i).sys.mode(1).D(:,m_u+m_d+m_w+1:m_u+m_v+m_w+m_d)), -kron(eye((T)), modes(j).sys.mode(1).D(:,m_u+m_d+m_w+1:m_u+m_v+m_w+m_d)); ...
                  -kron(eye((T)), modes(i).sys.mode(1).D(:,m_u+m_d+m_w+1:m_u+m_v+m_w+m_d)), kron(eye((T)), modes(j).sys.mode(1).D(:,m_u+m_d+m_w+1:m_u+m_v+m_w+m_d))];
    end
end

% Q denotes No. of pairs
Q=factorial(N)/(factorial(2)*factorial(N-2));

% Concatenate matrices further to obtain a concise optimization problem
for q=1:Q
    H_x{q}=Pbar_x*[M_x{q} Gamma_xd{q} Gamma_xw{q} zeros(size(Gamma_xw{q},1),size(Gamma_xw{q},2)/m_w*m_v)];
    H_y{q}=Pbar_y{q}*[M_y{q} Gamma_yd{q} Gamma_yw{q} zeros(size(Gamma_yw{q},1),size(Gamma_xw{q},2)/m_w*m_v)];
    
    h_x{q}=pbar_x - Pbar_x*ftilde_x{q}-Pbar_x*Gamma_xu{q}*u(:);
    h_y{q}=pbar_y{q} - Pbar_y{q}*ftilde_y{q}-Pbar_y{q}*Gamma_yu{q}*u(:);
    
    % Matrix for lumped uncertainty
    H_xbar{q} = blkdiag(kron(eye(2),bounds.P_0),Qbar_d{q},Qbar_w{q},Qbar_v{q});
    h_xbar{q} = [kron(ones(2,1),bounds.p_0); qbar_d{q};qbar_w{q};qbar_v{q}];
    
    % Matrix for separability condition
    G{q} = E{q}*[Abar_conc{q} Gamma_d_conc{q} Gamma_w_conc{q} zeros(size(Gamma_d_conc{q},1),size(Gamma_w_conc{q},2)/m_w*m_v)]...
           +[zeros(size(F_v{q},1),size(Abar_conc{q},2)),F_d{q},zeros(size(F_v{q},1),size(F_v{q},2)/m_v*m_w),F_v{q}];
    
    % Final matrices for optimizaiton problem
    % Inner problem: 
    % --------------------------
    % min   delta 
    % s.t.  R*x_bar <= [0,1]'*delta+r-S*u
    %       Phi*x_bar <= phi
    % --------------------------
    
    % Matrix for constraint being implicitly dependent on u
    Phi{q} = [H_xbar{q}];
    phi{q} = [h_xbar{q}];
 
    % Matrix for constraint being explicitly dependent on u
    % Without ego car's responsibility
    R{q} = [H_y{q}; G{q}];
    r{q}=[pbar_y{q}-Pbar_y{q}*ftilde_y{q}; -E{q}*ftilde{q}];
    S{q}=[Pbar_y{q}*Gamma_yu{q}; E{q}*Gamma_u_conc{q}+F_u{q}];
   
    % Outer problem (responsibility of controlled input):
    % --------------------------
    % R_x*x_bar <= r_x-S_x*u, for all x_abr: Phi_x*x_bar <= phi_x
    % --------------------------
    
    % Matrices for responsibility of controlled input
    R_x{q} = H_x{q};
    r_x{q} = pbar_x - Pbar_x*ftilde_x{q};
    S_x{q} = Pbar_x*Gamma_xu{q};  
    
    % Uncertainty
    Phi_x{q} = [H_xbar{q}; H_y{q}];
    % Pbar_y{q}*Gamma_yu{q} = 0 here (Assumption 1 is satisfied)
    phi_x{q} = [h_xbar{q}; pbar_y{q}-Pbar_y{q}*ftilde_y{q}-Pbar_y{q}*Gamma_yu{q}*u(:)]; 
     
end


%% To check whether the set of uncertaities is not empty

% Construct constraints in optimizaiton problem
constr = [];
constr = [constr;Qbar_u*u(:)<=qbar_u ];  % constraints on inputs

Q = factorial(N)/(factorial(2)*factorial(N-2)); % No. of mode pair

rho = 2*p*(T); % Total No. of output after concatenation

for q = 1:Q
    sigma = size(R{q},1);
    kappa = size(phi{q},1);
    eta = size(Phi{q},2);
 
    % Without the constraints for our/ego car's responsibility in the inner optimization problem
    xi=size([H_y{q}],1);
    
    % Declare decision variables
    mu1{q} = sdpvar(kappa,1,'full');
    mu2{q} = sdpvar(xi,1,'full');
    g1{q} = sdpvar(kappa,1,'full');
    g2{q} = sdpvar(xi,1,'full');
    mu3{q} = sdpvar(rho,1,'full');
    g3{q} = sdpvar(rho,1,'full');
    xbar{q} = sdpvar(eta,1,'full');
    delta{q} = sdpvar(1,1,'full');
    
    %% Constraint of outer problem: responsibility of controlled input 
    beta = size(R_x{q},1);
    psi = size(phi_x{q},1);
    Pi{q} = sdpvar(psi,beta,'full'); % dual variable
    
    % Constraint of outer problem: responsibility of controlled input
    % Robust counterpart
    constr_respond_ctr = [Pi{q}'*phi_x{q} <= r_x{q}-S_x{q}*u(:); ...
                          Pi{q}'*Phi_x{q} == R_x{q}; ...
                          Pi{q} >= 0];             
    constr = [constr; constr_respond_ctr];
    
    %% Constraint of inner problem: separability condition
    
    % Constraints on separation threshold
    constr = [constr; delta{q} >= epsi]; % Separability condition
    
    % KKT condition for inner problem
    %
    % Primal
    constr = [constr; R{q}(1:xi,:)*xbar{q} <= r{q}(1:xi)-S{q}(1:xi,:)*u(:); Phi{q}*xbar{q}<=phi{q}];
    constr = [constr; R{q}(1+xi:xi+rho,:)*xbar{q} <= ones(rho,1)*delta{q}+r{q}(1+xi:xi+rho)-S{q}(1+xi:xi+rho,:)*u(:)];
    
    % Dual
    constr = [constr; mu1{q}(:) >= 0; mu2{q}(:) >= 0; mu3{q}(:) >= 0]; 
    
    % Stationarity
    for m=1:eta
        constr = [constr; -mu1{q}(:)'*Phi{q}(:,m) == mu2{q}(:)'*R{q}(1:xi,m)+ mu3{q}(:)'*R{q}(1+xi:xi+rho,m)];
    end
    constr = [constr; sum(mu3{q}(:)) == 1];

    % Complementary slackness
    %
    % ---- use Big-M to reformulate bilinear constraint----
%     b1{q} = binvar(kappa,1,'full');
%     b2{q} = binvar(xi,1,'full');
%     b3{q} = binvar(rho,1,'full');
%     MM = 50;
%     MM1 = 50;
%     
%     for j=1:kappa
%         constr = [constr; g1{q}(j)-phi{q}(j) == -Phi{q}(j,:)*xbar{q}]; % need to be modifed according to big-M and complemets
%         constr = [constr; complements(mu1{q}(j) >= 0, g1{q}(j) >= 0)];
%         constr = [constr; mu1{q}(j) <= b1{q}(j,1)*MM; -(1-b1{q}(j,1))*MM1 <= g1{q}(j)];
%     end
%     
%     for i=1:xi
%         constr = [constr; g2{q}(i)+R{q}(i,:)*xbar{q}-r{q}(i) == -S{q}(i,:)*u(:)]; % need to be modifed according to big-M and complemets
%         constr = [constr; complements(mu2{q}(i) >= 0,g2{q}(i) >= 0)];
%         constr = [constr; mu2{q}(i) <= b2{q}(i,1)*MM; -(1-b2{q}(i,1))*MM1 <= g2{q}(i)];
%     end
%     
%     for i=1:rho
%         constr = [constr; g3{q}(i)+R{q}(xi+i,:)*xbar{q} == delta{q}+r{q}(xi+i)-S{q}(xi+i,:)*u(:)]; % need to be modifed according to big-M and complemets
%         constr = [constr; complements(mu3{q}(i) >= 0, g3{q}(i) >= 0)];
%         constr = [constr; mu3{q}(i) <= b3{q}(i,1); -(1-b3{q}(i,1))*MM1 <= g3{q}(i)];
%     end
    
    % ---- use SOS-1 to reformulate bilinear constraint----
    for j=1:kappa
        constr = [constr; g1{q}(j)+phi{q}(j) == Phi{q}(j,:)*xbar{q}; sos1([mu1{q}(j),g1{q}(j)],[1,10])];
    end
    
    for i=1:xi
        constr = [constr; g2{q}(i)-R{q}(i,:)*xbar{q}+r{q}(i)==S{q}(i,:)*u(:); sos1([mu2{q}(i),g2{q}(i)],[1,10])];
    end
    
    for i=1:rho
        constr = [ constr; g3{q}(i)-R{q}(xi+i,:)*xbar{q} == -delta{q}-r{q}(xi+i)+S{q}(xi+i,:)*u(:); sos1([mu3{q}(i),g3{q}(i)],[1,10]) ];
        constr = [constr; mu3{q}(i) <= 1];
    end
end

% Optimization
% Set objective function based on different norm
if (NNorm == 1 || NNorm == inf)
	obj = norm(u(:),NNorm);
elseif (NNorm == 3)
    obj = norm(u(:),1) + 2*norm(u(:),inf);
elseif NNorm == 2
	obj = norm(u(:),NNorm); % (u(:)'*u(:));  % 
end

options = sdpsettings('verbose', 1, 'debug', 1, 'solver', 'gurobi');
t = tic;
sol = optimize(constr,obj,options);
toc(t);
% Get optimal separating input
if sol.problem == 0
    ustar = value(u);
else
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
    ustar = value(u); % ustar = 0;
end
