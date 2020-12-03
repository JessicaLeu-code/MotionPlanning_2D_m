%% CHOMP generate functions 

clear all

% number of states
dim = 2;

% time steps
H = 5;

% sampling time
dt = 0.5;

%% %%%%%%%%%%%%%%  distance function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
x = sym('x',[dim 1]);

% Obstacles
c = sym('c',[dim 1]);

% margin 
syms m

% distance function
D = norm(x-c)-m;
%Df =  matlabFunction(D,'File','D_f');

% gradient 
dD = gradient(D,x)
%dDf =  matlabFunction(dD,'File','dD_f');


