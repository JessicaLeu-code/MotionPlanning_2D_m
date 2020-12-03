% generate functions 

clear all

% number of states
dim = 2;

% time steps
H = 5;

% sampling time
dt = 0.5;

% %% %%%%%%%%%%%%%%  distance function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % states
% x = sym('x',[dim 1]);
% 
% % Obstacles
% A = sym('A',[dim dim]);
% c = sym('c',[dim 1]);
% 
% % margin 
% syms m
% 
% % distance function
% g = -(x-c)'*A*(x-c)+m^2;
% %gf =  matlabFunction(g,'File','g_f');
% 
% % gradient 
% dg = gradient(g,x)
% %dgf =  matlabFunction(dg,'File','dg_f');


% %% %%%%%%%%%%%%%%  cost function   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% goal
goal = sym('g',[dim 1]);
% initial state
x0 = sym('x0',[dim 1]);

% inputs
u = sym('u',[dim*H 1]); 

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

% initialization
u0 = zeros(dim*H,1);
x_ = Aaug*x0 + Baug*u0;


%% cost function
Q = 10*eye(dim);                                % State cost
P = 100*eye(dim);                                % Terminal cost 
R = 50*eye(dim);                                % Input cost

Q_ = kron(eye(H-1),Q);
Q_ = blkdiag(Q_, P);
Raug = kron(eye(H),R);

% quadratic cost
Qaug = Baug'*Q_*Baug;
alpha = 1/max(svd(Qaug));

% linear cost 
g_aug = kron(ones(H,1),goal);
paug = [-(-x0'*Aaug'+g_aug')*Q_*Baug]';

% constant cost
caug = x0'*Aaug'*Q_*Aaug*x0;
cost = u'*Qaug*u + 2*paug'*u + caug;
%costf =  matlabFunction(cost,'File','cost2_f5');
% tic
% for i = 1 :100
     dcost = gradient(cost,u);
     subs(dcost, [x0; u], [[1 2]'; zeros(dim*H,1)]);
% end
% toc

% gradient 
dcost = gradient(cost,u);
%dcostf =  matlabFunction(dcost,'File','dcost2_f5');
