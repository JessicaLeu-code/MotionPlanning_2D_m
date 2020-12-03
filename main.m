%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSGCSF/CSF/CHOMP for xy-planning 
% Robot model: point mass (xy intergrater)
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all
Time_= [];
terminal_cost = [];
iter = [];
num_iter_chomp =0;
num_iter_cfs = 0;
num_iter_sqp = 0;
num_trials =1;

Time_rrtcfs = zeros(num_trials,1);
Time_rrt = zeros(num_trials,1);
num_iter =1;

%% Figures
figure(1); hold on
%% Environment
% Obstacle
obs{1}.shape = 'circle';
obs{1}.A = [1 0 ; 0 1];
obs{1}.c = [0.25, 0.02]';
obs{1}.m = 0.03999;
obs{1}.epsilon = 0.00;
obs{1}.v = [0;-0.0];

obs{2}.shape = 'circle';
obs{2}.A = [1 0 ; 0 1];
obs{2}.c = [0.06, 0.01]';
obs{2}.m = 0.02;
obs{2}.epsilon = 0.00;
obs{2}.v = [0;0];

obs{3}.shape = 'circle';
obs{3}.A = [1 0 ; 0 1];
obs{3}.c = [0.15, -0.03]';
obs{3}.m = 0.0399;
obs{3}.epsilon = 0.00;
obs{3}.v = [0;0];

% obs{4}.shape = 'rectangle';
% obs{4}.c = [0.2, 0.0]';
% obs{4}.w = 0.01;
% obs{4}.l = 0.05;
% obs{4}.poly = [obs{4}.c(1)-obs{4}.w/2  obs{4}.c(1)+obs{4}.w/2  obs{4}.c(1)+obs{4}.w/2  obs{4}.c(1)-obs{4}.w/2;
%                obs{4}.c(2)-obs{4}.l/2  obs{4}.c(2)-obs{4}.l/2  obs{4}.c(2)+obs{4}.l/2  obs{4}.c(2)+obs{4}.l/2];
% obs{4}.m = 0.02;
% obs{4}.epsilon = 0.02;
% obs{4}.v = [0;0];

num_obs = size(obs,2);

%% Robot
x0 = [0.4 0.0]';
% number of states
nstate = 2;
% number of inputs
dim = 2;

%% Initialization
H = 30;                                    % time steps
dt = 0.5;                                  % sampling time

%% Goal
goal = [0;0];

%% View and Plot
view_area = [-0.1 0.5 -0.2 0.2 0 0.1];
figure(1)
[plot_marker] = plot_obs_2D(obs,view_area);

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
R = 50*eye(dim);                                % Input cost

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
lambda = 500;
sys_info.lambda = lambda;
% solver setups
sys_info.options = optimoptions('quadprog','Display','off');
sys_info.epsilon_O = 1e-3;
sys_info.MAX_O_ITER = 200;
sys_info.float = 65.0783;
% xref
aa = linspace(x0(1),goal(1),H+1);
bb = linspace(x0(2),goal(2),H+1);
xref = [];
for i = 1:H+1
    xref = [xref; aa(i);bb(i)];
end

    
%% Solve
% baseline
eval = EVAL(sys_info);
[Cost_b,u_b] = eval.get_Cost_b();

pt = [];                                        % plot results
% CHOMP
tic
self = CHOMP(obs,sys_info);
self = self.optimizer();
% self = STOMP(obs,sys_info);
% self = self.optimizer();
Time_(1) = toc; 
self.x_ = [x0;self.x_];
Cost{1} = self.eval.cost_all;
terminal_cost(1) = self.eval.cost_all(end);
e_u_all{1} = self.eval.e_u_all;
e_cost_all{1} = self.eval.e_cost_all;
iter(1) = self.iter_O;
% planning result
pt(1) = plot(self.x_(1:2:end),self.x_(2:2:end),'-r*');
cost_chomp(1) = norm(self.x_(5:end)-g_aug)/H;
num_iter_chomp = num_iter_chomp+self.iter_O;


% PSGCFS
tic
self = PSGCFS(obs,sys_info);
self = self.optimizer();
Time_(2)= toc ;
Cost{2} = self.eval.cost_all;
terminal_cost(2) = self.eval.cost_all(end);
e_u_all{2} = self.eval.e_u_all;
e_cost_all{2} = self.eval.e_cost_all;
iter(2) = self.iter_O;
% planning result
pt(2) = plot(self.x_(1:2:end),self.x_(2:2:end),'g*');

% CFS
sys_info.solver = 'CFS';
tic
self = CFS(obs,sys_info);
self = self.optimizer();
Time_(3)= toc ;
Cost{3} = self.eval.cost_all;
terminal_cost(3) = self.eval.cost_all(end);
e_u_all{3} = self.eval.e_u_all;
e_cost_all{3} = self.eval.e_cost_all;
iter(3) = self.iter_O;
% planning result
num_iter_cfs = num_iter_cfs+self.iter_O;
cost_cfs(3) = norm(self.x_(3:end)-g_aug)/H;

% planning result
pt(3) = plot(self.x_(1:2:end),self.x_(2:2:end),'b*');
pt(4) = plot(goal(1),goal(2),'k*');
xlabel('x [m]')
ylabel('y [m]')
legend([pt plot_marker(1)],'CHOMP','PSGCFS','CFS','Goal','Obstacle','Location','northoutside','Orientation','horizontal')
axis(view_area)
axis equal
view(2)


%% Plot 
% cost calculation result
figure
plot(Cost{1},'r','LineWidth',1.5)
hold on 
plot(Cost{2},'g','LineWidth',1.5)
plot(Cost{3},'b','LineWidth',1.5)
yline(Cost_b,'-k');
yline(Cost{1}(end),'-r');
yline(Cost{2}(end),'-g');
yline(Cost{3}(end),'-b');
legend('CHOMP','PSGCFS','CFS')
xlabel('Iteration')
ylabel('Cost')
