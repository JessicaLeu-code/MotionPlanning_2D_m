%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFS for xy-planning 
% Robot model: xy intergrater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

fighandle = [];
fighandle(1) = figure(1); hold on;

N = 21;
path = [linspace(-2,1,N);
        linspace(0,0,N)];
    
dt = 1;nobj = 1;
obs = {};
obs{1}.v = [0;0];
dx = -1;
dy = 0;
obs{1}.poly = [0.2+dx 0.8+dx 0.81+dx 0.19+dx; -0.30001+dy -0.30001+dy 1+dy 1+dy];

nstep = size(path,2);
dim = 2;
MAX_iter = 10;

%% The cost function
% The distance metric between the original path and the new path
Q1 = zeros(nstep*dim);
Q1((nstep-1)*dim+1:end,(nstep-1)*dim+1:end) =  eye(dim)*1000;
Q1(1:dim,1:dim) =  eye(dim)*1000;
% The velocity
Vdiff = eye(nstep*dim)-diag(ones(1,(nstep-1)*dim),dim);
Q2 = Vdiff(1:(nstep-1)*dim,:)'*Q1(1+dim:end,1+dim:end)*Vdiff(1:(nstep-1)*dim,:);
% The accelaration
Adiff = Vdiff-diag(ones(1,(nstep-1)*dim),dim)+diag(ones(1,(nstep-2)*dim),dim*2);
Q3 = Adiff(1:(nstep-2)*dim,:)'*Adiff(1:(nstep-2)*dim,:);
% The weight
c = [0,2,2];
% The total costj
Qref = 1*(Q1*c(1)+Q2*c(2)+Q3*c(3));
Qabs = 0*Q3*c(3);
%% Extended cost
Mcurv = eye(nstep);
Mcurv(nstep,nstep) = 5;
Vcurv = eye(nstep)-diag(ones(1,nstep-1),1);
Acurv = Vcurv-diag(ones(1,(nstep-1)),1)+diag(ones(1,(nstep-2)),2);
Qcurv = 5*Mcurv;%+Vcurv(1:nstep-1,:)'*Vcurv(1:nstep-1,:)+Acurv(1:(nstep-2),:)'*Acurv(1:(nstep-2),:);
%% The boundary constraint
Aeq = zeros(4*dim,nstep*dim+nstep);
Aeq(0*dim+1:1*dim,1:dim) = eye(dim);
Aeq(1*dim+1:2*dim,(nstep-1)*dim+1:nstep*dim) = eye(dim);
Aeq(2*dim+1:3*dim,1:2*dim) = [-eye(dim) eye(dim)];
Aeq(3*dim+1:4*dim,(nstep-2)*dim+1:nstep*dim) = [-eye(dim) eye(dim)];
beq = [path(:,1);path(:,end);path(:,2)-path(:,1);path(:,end)-path(:,end-1)];

refpath = [];
for i=1:nstep
    refpath = [refpath;path(:,i)];
end
oripath = refpath;
refinput = ones(1,nstep);

%% The Iteration
tic
for k = 1:MAX_iter
FEAS = 1;
%% The constraint
Lstack = []; Sstack = []; margin = 0.2;
for i=1:nstep
    for j=1:nobj
        poly = obs{j}.poly+obs{j}.v*ones(1,4)*dt*i;
        [L,S,d] = d2poly(refpath((i-1)*dim+1:i*dim)',poly');
        Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,nstep)];
        Sstack = [Sstack;S-margin];
    end

end

%% QP
Qe = blkdiag(Qref+Qabs,0*Qcurv);

soln = quadprog(Qe,[-Qref*oripath;zeros(nstep,1)],Lstack,Sstack,Aeq,beq);
pathnew = soln(1:dim*nstep);
refinput = soln(dim*nstep+1:end);


figure(fighandle(1));
plot(pathnew(1:dim:end),pathnew(2:dim:end),'-*','color',[1-k/MAX_iter,1-k/MAX_iter,1-k/MAX_iter])
if norm(refpath-pathnew) < 0.1
    disp(['converged at step ',num2str(k)]);
    break
end
refpath = pathnew;
end
%%
time = toc
disp(['final cost: ']);
cost_curv(pathnew,oripath,Qref,Qabs,Qcurv,nstep)

%%
figure(fighandle(1));
plot(pathnew(1:dim:end),pathnew(2:dim:end),'k')
ob = Polyhedron('V',obs{1}.poly');
ob.plot('color','g');
axis equal
grid off
box on
legend('Iter1','Iter2','Iter3','Iter4','Iter5')