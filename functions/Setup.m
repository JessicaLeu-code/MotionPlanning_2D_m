function  [sys_info_out,self_out] = Setup(dt,x0,goal,obs,SOLVER,self)

persistent N sys_info Aaug Baug Q_ H self_

if isempty(N)
    % initialize all necessary parameters     
    N       = uint16(1);
    
    nstate = 2;                                % number of states
    dim = 2;                                   % number of inputs
    H = 30;                                    % time steps
                                     
    
    %% dynamics
    Ad = eye(dim);
    Bd = eye(dim)*dt;
    Aaug = [];
    Baug = zeros(H*dim,H*dim);
    for i=1:H
        Aaug = [Aaug;Ad^i];
        for j=1:i
            Baug((i-1)*dim+1:i*dim,(j-1)*dim+1:j*dim)=Ad^(i-j)*Bd;
        end    
    end

    %% cost function
    Q = 10*eye(dim);                                % State cost
    P = 1000*eye(dim);                              % Terminal cost 
    R = 10*eye(dim);                                % Input cost

    Q_ = kron(eye(H-1),Q);
    Q_ = blkdiag(Q_, P);
    Raug = kron(eye(H),R);

    % quadratic cost
    Qaug = Baug'*Q_*Baug+10*Raug;
    alpha = 1/max(svd(Qaug));

    % linear cosr 
    g_aug = kron(ones(H,1),goal);
    paug = [-(-x0'*Aaug'+g_aug')*Q_*Baug]';

    % constant cost
    caug = (-x0'*Aaug'+g_aug')*Q_*(-x0'*Aaug'+g_aug')';

    %% System info
    sys_info.num_obs = max(size(obs));
    sys_info.H = H;
    sys_info.dim = dim;
    sys_info.dt = dt;
    sys_info.nstate = nstate;
    sys_info.x0 = x0;
    sys_info.goal = goal;
    sys_info.Aaug = Aaug;
    sys_info.Baug = Baug;
    sys_info.Qaug = Qaug;
    sys_info.paug = paug;
    sys_info.caug = caug;

    % limits
    sys_info.in_max = [0.2 0.2]'*dt;                % velocity limits
    sys_info.alpha = alpha;
    % CHOMP penalty ratio
    lambda = 100;
    sys_info.lambda = lambda;
    % solver setups
    sys_info.options = optimoptions('quadprog','Display','off');
    sys_info.epsilon_O = 1e-3;
    sys_info.MAX_O_ITER = 100;
    
    
    % Initialize solvers
    switch(SOLVER)
        case 'CHOMP'
            self_ = CHOMP(obs,sys_info);
        case 'PSGCFS'
            self_ = PSGCFS(obs,sys_info);
        case 'CFS'
            sys_info.solver = SOLVER;
            self_ = CFS(obs,sys_info);
    end
else
    self_ = self;
    self_.iter_O = 1;    
    g_aug = kron(ones(H,1),goal);
    paug = [-(-x0'*Aaug'+g_aug')*Q_*Baug]';
    % constant cost
    caug = (-x0'*Aaug'+g_aug')*Q_*(-x0'*Aaug'+g_aug')';
    sys_info.x0 = x0;
    sys_info.paug = paug;
    sys_info.caug = caug;
    sys_info.num_obs = max(size(obs));
    self_.sys_info = sys_info;
    self_.obs = obs;
end

sys_info_out = sys_info;
self_out = self_;

end