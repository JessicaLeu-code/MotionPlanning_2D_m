%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSGCSF/CSF for xy-planning
% (MPC implementation)
% Robot model: point mass (xy intergrater)
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
Time_= [];
terminal_cost = [];
iter = [];

%% Solver
SOLVER ='PSGCFS';  % will reach the goal
%SOLVER = 'CFS';   % will stuck infornt of the rectangular obstacle

%% Environment
% Obstacle
ran_=0;

obs{1}.shape = 'circle';
obs{1}.A = [1 0 ; 0 1];
obs{1}.c = [0.5, -0.0]'+ran_*0.005*(rand(2,1)-0.5);
obs{1}.m = 0.0399;%+ran_*0.01*(rand(1)-0.5);
obs{1}.epsilon = 0.00;
obs{1}.v = [0;-0.0];

obs{2}.shape = 'rectangle';
obs{2}.c = [0.15, 0]';
obs{2}.w = 0.04;
obs{2}.l = 0.05;
obs{2}.poly = [obs{2}.c(1)-obs{2}.w/2  obs{2}.c(1)+obs{2}.w/2  obs{2}.c(1)+obs{2}.w/2  obs{2}.c(1)-obs{2}.w/2;
               obs{2}.c(2)-obs{2}.l/2  obs{2}.c(2)-obs{2}.l/2  obs{2}.c(2)+obs{2}.l/2  obs{2}.c(2)+obs{2}.l/2];
obs{2}.m = 0.02;
obs{2}.epsilon = 0.02;
obs{2}.v = [0.01;0];


%% Robot
x0 = [0.6 0.0]';
% number of states
nstate = 2;
% number of inputs
dim = 2;

%% Initialization
H = 30;                                    % time steps
dt = 0.2;                                  % sampling time

%% Goal
goal = [0;0];

%% View and Plot
view_area = [-0.1 0.9 -0.3 0.3 0 0.1];
figure(1)
[plot_marker] = plot_obs_2D(obs,view_area);
start = plot(x0(1),x0(2),'b*');

%% MPC loop start
self=[];
steps = 150;
ss = 1;
while ss<80 && x0(1)>0.01
%% Setup sys_info
[sys_info,self] = Setup(dt,x0,goal,obs,SOLVER,self);
    
%% solve
pt = [];                                    % plot results
tic
self = self.optimizer();
Time_(ss)= toc;                       
Cost{ss} = self.eval.cost_all;
terminal_cost(ss) = self.eval.cost_all(end);
e_u_all{3,ss} = self.eval.e_u_all;
e_cost_all{ss} = self.eval.e_cost_all;
iter(ss) = self.iter_O;
% planning result
pt(3) = plot(self.x_(1:2:end),self.x_(2:2:end),'b*');
pt(4) = plot(goal(1),goal(2),'c*');
xlabel('x [m]')
ylabel('y [m]')
axis(view_area)
axis equal
view(2)
pause(0.1)
delete(pt)

%% Next step
x0 = self.x_(self.sys_info.dim+1:self.sys_info.dim*2);
ref_length = 29;
self.x_ = [self.x_(self.sys_info.dim+1:self.sys_info.dim*ref_length);kron(ones(self.sys_info.H-(ref_length-1)+1,1),self.x_(self.sys_info.dim*ref_length+1:self.sys_info.dim*(ref_length+1)))];
plot(x0(1),x0(2),'ro')

[obs] = update_obs(obs,sys_info.dt);
delete(plot_marker);
[plot_marker] = plot_obs_2D(obs,view_area);

ss = ss + 1;

end

