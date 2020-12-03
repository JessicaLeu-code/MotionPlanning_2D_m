%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOMP class file for xy-planning 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef STOMP
   properties
       % setup parameters
       obs cell
       sys_info struct
       nn {mustBeNumeric}
       SOLVER = 'CHOMP'                    % 'CHOMP' or 'CHOMP-HMC'
             
       % outputs
       u {mustBeNumeric}                   % velocity
       x_ {mustBeNumeric}                  % trajectory
       
       % solver settings
       num_samples = 20;
       A {mustBeNumeric}                   % Qaug
       R_inv {mustBeNumeric}               % inv(A'*A)
       h = 10;
       
                   
       % evaluation setting
       eval EVAL
       
       % calculation statistics       
       iter_O = 1;
       total_iter = 0;
       
   end
   
   methods
       % _init_
       function self = STOMP(val,val2,varargin)
           % get problem info
           self.obs = val;
           self.sys_info = val2;                
           self.nn = self.sys_info.H*self.sys_info.dim;
           if ~isempty(varargin)
                    self.SOLVER = varargin{1};                    
           end 
           % initialize solution
           self.u = zeros(self.nn,1);           
           self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
           self.A = self.sys_info.Qaug;
           self.R_inv = inv(self.A'*self.A);
           % initialize evaluation functions
           self.eval = EVAL(val2);
       end
       
       % main function
       function self = optimizer(self)   
           self.eval.cost_new = self.getCost();
           while self.eval.stop_outer(self.iter_O)~=true
               self.eval.cost_old = self.eval.cost_new;
               self.eval.u_old = self.u;
               self.u = self.u - self.get_du();
               self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
               self.eval.cost_new = self.eval.get_cost(self.u);%self.getCost();
               self.eval = self.eval.store_result(self.u);
               self.iter_O = self.iter_O+1;
           end               
                  
       end
       
       %%%%%%%%%%%%%%%% supplementary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function  [dU] = dU_f(self)
           self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
           dco = zeros(self.sys_info.H*self.sys_info.dim,1);
           for j= 1:self.sys_info.num_obs
               for i = 1:self.sys_info.H
                   dco = dco + [dc_f(self.x_((i-1)*self.sys_info.dim+1:i*self.sys_info.dim),self.obs{j}.c,self.obs{j}.m,self.obs{j}.epsilon)'*self.sys_info.Baug((i-1)*self.sys_info.dim+1:i*self.sys_info.dim,:)]';  
               end
           end
           dU = (2*(self.sys_info.Qaug*self.u+self.sys_info.paug)+self.sys_info.lambda*dco);
       end
       
       function Ux = getCost(self)           
           cobs = 0;
           for j= 1:self.sys_info.num_obs
               for ii = 1:self.sys_info.H
                   cobs = cobs + c_f(self.x_((ii-1)*self.sys_info.dim+1:ii*self.sys_info.dim),self.obs{j}.c,self.obs{j}.m,self.obs{j}.epsilon);
               end
           end
           Ux = self.eval.get_cost(self.u)+ self.sys_info.lambda*cobs;
       end
       
              
       function [du] = get_du(self)
           % initialize
           epsilon = zeros(self.nn,self.num_samples);
           u_per = zeros(self.nn,self.num_samples);
           x_per = zeros(self.nn,self.num_samples);
           Sth = zeros(self.sys_info.H,self.num_samples);
           P = zeros(self.sys_info.H,self.num_samples);
           % get samples
           for k = 1:self.num_samples
               epsilon(:,k) = mvnrnd(zeros(self.nn,1),self.R_inv);
               u_per(:,k) = self.u + epsilon(:,k);
               x_per(:,k) = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*u_per(:,k);
           end
           % get weight
           for i = 1:self.sys_info.H
               for j= 1:self.sys_info.num_obs
                   for k = 1:self.num_samples
                       Sth(k,i) = self.sys_info.paug((i-1)*self.sys_info.dim+1:i*self.sys_info.dim)'*u_per((i-1)*self.sys_info.dim+1:i*self.sys_info.dim,k)... 
                           + c_f(x_per((i-1)*self.sys_info.dim+1:i*self.sys_info.dim,k),self.obs{j}.c,self.obs{j}.m,self.obs{j}.epsilon);
                   end
               end
               Smax = max(Sth(:,i));
               Smin = min(Sth(:,i));
               for k = 1:self.num_samples
                   P(i,k) = exp(-self.h*(Sth(i,k)-Smin)/(Smax-Smin));                   
               end
               P(i,:) = P(i,:)./sum(P(i,:));
           end
           
           du_ = zeros(self.sys_info.H,1);
           for i = 1:self.sys_info.H               
                   for k = 1:self.num_samples
                       du_((i-1)*self.sys_info.dim+1:i*self.sys_info.dim) = P(i,k)*epsilon((i-1)*self.sys_info.dim+1:i*self.sys_info.dim,k);
                   end
           end
           du = (self.R_inv./(-min(self.R_inv)))*du_; 
       end
             
       %%%%%%%%%%%%%%% HMC functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   end
end