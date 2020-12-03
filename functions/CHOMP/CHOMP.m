%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOMP class file for xy-planning 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef CHOMP
   properties
       % setup parameters
       obs cell
       sys_info struct
       nn {mustBeNumeric}
       SOLVER = 'CHOMP'                    % 'CHOMP' or 'CHOMP-HMC'
             
       % outputs
       u {mustBeNumeric}                   % velocity
       x_ {mustBeNumeric}                  % trajectory
                   
       % evaluation setting
       eval EVAL
       
       % calculation statistics       
       iter_O = 1;
       total_iter = 0;
       
   end
   
   methods
       % _init_
       function self = CHOMP(val,val2,varargin)
           % get problem info
           self.obs = val;
           self.sys_info = val2;                
           self.nn = self.sys_info.H*self.sys_info.dim;
           % initialize solution
           self.u = zeros(self.nn,1);
           if ~isempty(varargin)
               self.SOLVER = varargin{1};
               self.x_ = varargin{2};
           else
               self = self.get_path();
           end        
           % initialize evaluation functions
           self.eval = EVAL(val2);
       end
       
       % main function
       function self = optimizer(self)   
           self.eval.cost_new = self.getCost();
           switch self.SOLVER
               case 'CHOMP' 
                   while self.eval.stop_outer(self.iter_O)~=true
                       self.eval.cost_old = self.eval.cost_new;
                       self.eval.u_old = self.u;
                       self.u = self.u - (self.sys_info.alpha/self.iter_O)*self.dU_f();
                       self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
                       self.eval.cost_new = self.eval.get_cost(self.u);%self.getCost();
                       self.eval = self.eval.store_result(self.u);
                       self = self.get_path();
                       self.iter_O = self.iter_O+1;
                   end
               case 'CHOMP-HMC'
                   self = self.HMC();
                   
           end
                   
       end
       
       %%%%%%%%%%%%%%%% supplementary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function  [dU] = dU_f(self)
           self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
           dco = zeros(self.sys_info.H*self.sys_info.dim,1);
           for j= 1:self.sys_info.num_obs
               for i = 1:self.sys_info.H
                   if ~isempty(self.obs{j})
                       dco = dco + [dc_f(self.x_((i-1)*self.sys_info.dim+1:i*self.sys_info.dim),self.obs{j}.c,self.obs{j}.m,self.obs{j}.epsilon)'*self.sys_info.Baug((i-1)*self.sys_info.dim+1:i*self.sys_info.dim,:)]';  
                   end
               end
           end
           dU = (2*(self.sys_info.Qaug*self.u+self.sys_info.paug)+self.sys_info.lambda*dco);
       end
       
       function self = get_path(self)
           x0_ = self.sys_info.x0;
           self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
           self.x_ = [x0_;self.x_]; 
                      
       end
       
       function Ux = getCost(self)           
           cobs = 0;
           for j= 1:self.sys_info.num_obs
               for ii = 1:self.sys_info.H
                   if ~isempty(self.obs{j})
                        cobs = cobs + c_f(self.x_((ii-1)*self.sys_info.dim+1:ii*self.sys_info.dim),self.obs{j}.c,self.obs{j}.m,self.obs{j}.epsilon);
                   end
               end
           end
           Ux = self.eval.get_cost(self.u)+ self.sys_info.lambda*cobs;
       end
             
       %%%%%%%%%%%%%%% HMC functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function self = HMC(self)
           ee = self.sys_info.alpha;
           for k = 1:1
    
            % resample
            ut = self.u;
            %rt = (rand(dim*H,1)-0.5)*(alpha/k);
            rt = normrnd(0,1,[self.nn,1]);

            for kk = 1:600
                % leapfrog
                L = 15;                
                rt_ = rt - (ee/2)*self.dU_f(); 
                old_u = self.u;
                
                for ll = 1:L
                    self.u = ut + ee*rt_;
                    rt_ = rt_ -(ee)*self.dU_f();
                end
                rt_ee = rt_ -(ee/2)*self.dU_f();
                rt_ee = -rt_ee;               
                
                if self.P_f(rt,old_u,k)<self.P_f(rt_ee,self.u,k)
                    rt = rt_ee;
                    ut = self.u;                    
                else
                    pp = self.P_f(rt_ee,self.u,k)/self.P_f(rt,old_u,k);
                    
                    if rand<pp
                        rt = rt_ee;
                        ut = self.u;
                    else
                        self.u = old_u;
        %                 rt = (rand(dim*H,1)-0.5)*(alpha);
        %                 rt = norm(old_rt)*rt/norm(rt);
                        rt = normrnd(0,1,[self.nn,1]);
                    end
                end
                self.u_all = [ self.u_all self.u];
                self.eval.cost_all = [self.eval.cost_all self.getCost()];
            end           

           end
           self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*self.u;
           
       end
       
       function [Pfx,Hx] = P_f(self,rt,u,k)
           self.x_ = self.sys_info.Aaug*self.sys_info.x0 + self.sys_info.Baug*u;
           cobs = 0;
           for j= 1:self.sys_info.num_obs
               for ii = 1:self.sys_info.H
                   cobs = cobs + c_f(self.x_((ii-1)*self.sys_info.dim+1:ii*self.sys_info.dim),self.obs{j}.c,self.obs{j}.m,self.obs{j}.epsilon);
               end
           end
           Ux = 0.5*u'*self.sys_info.Qaug*u +self.sys_info.paug'*u + self.sys_info.caug + cobs;
           Kx = 0.5*rt'*rt;
           Hx = Ux+Kx;
           Pfx = exp(-(k/2000)*(Ux+Kx));
       end
       
   end
end