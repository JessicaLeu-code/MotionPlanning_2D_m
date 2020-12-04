%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RRT-CFS for xy-planning 
% Robot model: point mass (xy intergrater)
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
num_tri = 1;
Time_ = zeros(num_tri,1);
Time_rrt = zeros(num_tri,1);
num_iter =100;

for tri =1:num_tri
%% Environment
% Obstacle
ran_ = 1;                % add randomness

% % obstacle avoidance
% obs{1}.shape = 'circle';
% obs{1}.A = [1 0 ; 0 1];
% obs{1}.c = [0.5, -0.0]'+ran_*0.005*(rand(2,1)-0.5);
% obs{1}.m = 0.0399;%+ran_*0.01*(rand(1)-0.5);
% obs{1}.epsilon = 0.00;
% obs{1}.v = [0;-0.0];
% 
% obs{3}.shape = 'circle';
% obs{3}.A = [1 0 ; 0 1];
% obs{3}.c = [0.2, 0.02]'+ran_*0.1*(rand(2,1)-0.5);
% obs{3}.m = 0.0399+ran_*0.05*(rand(1)-0.5);
% obs{3}.epsilon = 0.00;
% obs{3}.v = [0;0];
% 
% obs{2}.shape = 'circle';
% obs{2}.A = [1 0 ; 0 1];
% obs{2}.c = [0.06, 0.01]'+ran_*0.01*(rand(2,1)-0.5);
% obs{2}.m = 0.02+ran_*0.01*(rand(1)-0.5);
% obs{2}.epsilon = 0.00;
% obs{2}.v = [0;0];

% Narrow passage problem
obs{1}.shape = 'rectangle';
obs{1}.c = [0.2, 0.439+0.115]';
obs{1}.w = 0.2;
obs{1}.l = 0.88;
obs{1}.poly = [obs{1}.c(1)-obs{1}.w/2  obs{1}.c(1)+obs{1}.w/2  obs{1}.c(1)+obs{1}.w/2  obs{1}.c(1)-obs{1}.w/2;
               obs{1}.c(2)-obs{1}.l/2  obs{1}.c(2)-obs{1}.l/2  obs{1}.c(2)+obs{1}.l/2  obs{1}.c(2)+obs{1}.l/2];
obs{1}.m = 0.01;
obs{1}.epsilon = 0.02;
obs{1}.v = [0;0];

obs{2}.shape = 'rectangle';
obs{2}.c = [0.2, -0.439+0.105]';
obs{2}.w = 0.2;
obs{2}.l = 0.88;
cc = obs{2}.c;
ww = obs{2}.w;
ll = obs{2}.l;
obs{2}.poly = [cc(1)-ww/2  cc(1)+ww/2  cc(1)+ww/2  cc(1)-ww/2;
               cc(2)-ll/2  cc(2)-ll/2  cc(2)+ll/2  cc(2)+ll/2];
obs{2}.m = 0.01;
obs{2}.epsilon = 0.02;
obs{2}.v = [0;0];

sys_info.num_obs = size(obs,2);

%% View and Plot
view_area = [-0.1 0.9 -0.3 0.3 0 0.1];
[plot_marker] = plot_obs_2D(obs,view_area);

%% Robot
% point mass
% number of states
nstate = 2;
% number of inputs
dim = 2;

%% Initialization
% initial position
x0 = [0.7; 0.0];
start = plot(x0(1),x0(2),'b*');

%% Goal
goalxyth = [0;0];
goal_th = goalxyth;
region_g = [0.02 0.02;]';

%% State constraints
region_s = [0.5;0.5];
sample_off = [0;0];

%% System info
% solver setups
sys_info.dim = dim;
sys_info.nstate = nstate;
sys_info.x0 = x0;
sys_info.goal_th = goal_th;


%% solver RRT
tic
% single-thread
% self = RRT_PM(obs,sys_info,region_g,region_s,sample_off,'RRT');
% self = self.RRTfind();

% multi-thread
self_ = {};
num_seed = 6;                      % depends on #cores of the computer 
path_fail = true(num_seed,1);
while all(path_fail)
    routeL = 1000*ones(num_seed,1);
    parfor i=1:num_seed 
    self_{i} = RRT_PM(obs,sys_info,region_g,region_s,sample_off,'RRT');
    self_{i} = self_{i}.RRTfind();
    path_fail(i) = self_{i}.fail;
        if self_{i}.fail ~=true
            routeL(i) = size(self_{i}.route,2);    
        end
    end
    num_iter = num_iter+1;
end
%%
[path_length,id] = min(routeL);
self = self_{id};
Time_rrt(tri)= toc;

%%%%%%%%%%%%% CFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%H = 30;
H = 30;%size(sampled_route,2)-1;           % time steps
dt = 0.5;                                  % sampling time
u0 = zeros(dim*H,1);                       % input initialization

horizon = H ;
wpTimes = (0:size(self.route,2)-1)*dt;
trajTimes = linspace(0,wpTimes(end),horizon+1);
sampled_route = cubicpolytraj(self.route,wpTimes,trajTimes);

%% Goal
goal = [0;0];

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
P = 5000*eye(dim);                              % Terminal cost 
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
% reference from RRT/RRT* result
xref = [];
for i = 1:H+1
    xref = [xref; sampled_route(:,i)];
end

% limits
sys_info.in_max = [0.2 0.2]'*dt;                % velocity limits
sys_info.alpha = alpha;
sys_info.solver = 'CFS';
sys_info.options = optimoptions('quadprog','Display','off');
% solver setups
sys_info.epsilon_O = 1e-3;
sys_info.MAX_O_ITER = 200;

%% CFS solver
self_cfs = CFS(obs,sys_info,xref);     
self_cfs = self_cfs.optimizer();
Time_(tri)= toc;
num_iter = num_iter+self_cfs.iter_O;
eval = EVAL(sys_info);
Cost_b(tri) = eval.get_Cost_b();



rrt_route= [];
for j =1:H
    rrt_route = [rrt_route; sampled_route(:,j+1) ];
end
cost_rrtcfs(tri) = norm(self_cfs.x_(3:end)-g_aug)/H;
cost_rrt(tri) = norm(rrt_route-g_aug)/H;
end

%% Plot
path_length = size(self.route,2);
pathimplemented = self.route(1:2,:);

pt(1) = plot(goalxyth(1),goalxyth(2),'*b')
pt(2) = plot(self.all_nodes(2,:),self.all_nodes(3,:),'ok');
route = plot(self.route(1,:),self.route(2,:),'*-r');

figure(1)
pt(3) = plot(self_cfs.x_(1:2:end),self_cfs.x_(2:2:end),'b*');
p_init = plot(x0(1),x0(2),'co');
p_goal = plot(goal_th(1),goal_th(2),'k*');
legend([pt p_goal p_init  plot_marker(1)],'RRT samples','RRT*','RRT*-CFS','Goal','Initial location','Obstacle','Location','northoutside','Orientation','horizontal','NumColumns',4)

axis(view_area)
axis([-0.1 0.9 -0.3 0.3])
axis equal
view(2)
