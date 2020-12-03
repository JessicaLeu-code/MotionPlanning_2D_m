%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSGCFS class file for xy-planning  
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PSGCFS
   properties
       % setup parameters
       obs cell
       sys_info struct
       nn {mustBeNumeric}       
      
       % outputs
       u {mustBeNumeric}                   % velocity
       x_ {mustBeNumeric}                  % trajectory
            
       % constraint
       Ainq {mustBeNumeric}
       binq {mustBeNumeric}       
       
       % inner sover settings
       MAX_I_ITER = 1;
       epsilon_I = 1e-3;
       num_sample = 10;
       u_batch {mustBeNumeric}
       
       % evaluation setting
       eval EVAL
       cost_old_I =1000;
       
       % calculation statistics
       iter_I = 1;
       iter_O = 1;
       total_iter = 0;
       
   end
   methods
       % _init_
       function self = PSGCFS(val,val2,varargin)
           % get problem info
           self.obs = val;
           self.sys_info = val2;                
           self.nn = self.sys_info.H*self.sys_info.dim;
           % initialize solution
           self.u = zeros(self.nn,1);
           if ~isempty(varargin)
               self.x_ = varargin{1};
           else
               self = self.get_path();
           end 
           
           self.u_batch = zeros(self.nn,self.MAX_I_ITER);
           % initialize evaluation functions
           self.eval = EVAL(val2);
                               
       end
       
       % main function
       function self = optimizer(self)           
           self.eval.cost_old = 10000;           
           self.eval.cost_new = self.eval.get_cost(self.u);
           while self.eval.stop_outer(self.iter_O)~=true                
               self.eval.cost_old = self.eval.cost_new;
               % get FS               
               %self.plot();
               self = self.get_Gamma_AB();
               % do PSG
               %%%%%%%%%%%% inner %%%%%%%%%%%%%
               self = self.inner_PSG();
               %self = self.inner_PSG_float();
%                self.eval.u_old = self.u;                    
%                self.eval.cost_old = self.eval.cost_new;
%                self = self.PSG_update();
%                % get cost
%                self.eval.cost_new = self.eval.get_cost(self.u);
%                % store results
%                self.eval = self.eval.store_result(self.u);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               self = self.get_path();
               % next step
               self.iter_O = self.iter_O+1; 
               % reset                             
               self.iter_I = 1;
               
           end           
       end  
              
       %%%%%%%%%%%%%%%% solver functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function self = inner_PSG(self)
           self.eval.u_old = self.u;
           while self.stop_inner()~=true
               self.eval.u_old = self.u;
               self.eval.cost_old = self.eval.cost_new;
               self = self.PSG_update();
               % get cost
               self.eval.cost_new = self.eval.get_cost(self.u);
           end
           self.eval = self.eval.store_result(self.u);
       end
        
       function self = inner_PSG_float(self)
           self.eval.u_old = self.u;
           while self.stop_inner()~=true
               %self.eval.cost_old = self.eval.cost_new;
               self = self.PSG_update_float();               
           end
           % avg from batch
           u_ = sum(self.u_batch,2)/size(self.u_batch,2);
           self = self.Projection(u_);
           % get cost
           self.eval.cost_new = self.eval.get_cost(self.u);
           %self.sys_info.float = self.eval.cost_new -0.1;
           self.eval = self.eval.store_result(self.u);
        end
       
       function self = PSG_update(self)
           u_ = self.u;
           %u_ = u_ - self.sys_info.alpha*(dcost2_f15(self.sys_info.goal,u_,self.sys_info.x0)+1*normrnd(0,0.1,[self.nn,1])./((self.iter_I)^5 + 1));
           randd = normrnd(0,pi*50,[self.nn,1])./((self.iter_O)^9 + 1);
           u_ = u_ - self.sys_info.alpha*(2*(self.sys_info.Qaug*self.u+self.sys_info.paug)+1*randd);%normrnd(0,0.01,[self.nn,1])./((self.iter_O)^5 + 1));
%            u_ = u_ - self.sys_info.alpha*(2*(self.sys_info.Qaug*self.u+self.sys_info.paug)+1*kron(ones(self.sys_info.H,1),normrnd(0,pi,[2,1]))./(self.iter_O^5 + 1));
%            self.u = u_;
%            self = self.get_path();
%            self = self.get_Gamma_AB();
           self = self.Projection(u_);           
           self.iter_I = self.iter_I +1;           
       end
       
       function self = PSG_update_float(self)
           u_ = self.u;
           randd = normrnd(0,pi*10000,[self.nn,1])./((self.iter_O)^5 + 1);
           u_ = u_ - self.sys_info.alpha*0.2*(2*(self.sys_info.Qaug*self.u+self.sys_info.paug)*(self.eval.cost_old-self.sys_info.float)+1*randd);%normrnd(0,0.01,[self.nn,1])./((self.iter_O)^5 + 1));
           self = self.Projection(u_);           
           self.iter_I = self.iter_I +1;
           self.u_batch(:,self.MAX_I_ITER) = self.u;
       end
       
%        function self = PSG_update2(self)
%            self.eval.u_old = self.u;
%            u_ = self.u;
%            temp_best = 1010101010;           
%            for i = 1:self.num_sample
%            randd = normrnd(0,pi*100,[self.nn,1])./((self.iter_O)^9 + 1);
%            u_ = u_ - self.sys_info.alpha*(2*(self.sys_info.Qaug*self.u+self.sys_info.paug)+1*randd);%normrnd(0,0.01,[self.nn,1])./((self.iter_O)^5 + 1));
%            H = eye(self.nn);
%            f = -u_;
%            [u_proj,fval,exitflag,output,lambda] = quadprog(H,f,self.Ainq,self.binq,[],[],[],[],[],self.sys_info.options);
%            temp_c = self.eval.get_cost(u_proj);
%                if temp_best > temp_c
%                    out_u = u_proj;
%                    temp_best = temp_c;
%                end
%            end
%            self.u = out_u;
%            %self = self.Projection(u_);           
%            self.iter_I = self.iter_I +1;           
%        end
              
       function self = Projection(self,u_)          
           H = eye(self.nn);
           f = -u_;
           [u_proj,fval,exitflag,output,lambda] = quadprog(H,f,self.Ainq,self.binq,[],[],[],[],[],self.sys_info.options);
           self.total_iter = self.total_iter + output.iterations;
           self.u  = u_proj; 
       end 
       
       function self = get_Gamma_AB(self)
           [Ainq_, binq_] = feasible_set(self.obs,self.x_,self.sys_info);
           self.Ainq = Ainq_;
           self.binq = binq_;
       end
       
       %%%%%%%%%%%%%%%% supporting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function self = get_path(self)
           x0_ = self.sys_info.x0;
           self.x_ = self.sys_info.Aaug*x0_ + self.sys_info.Baug*self.u;
           self.x_ = [x0_;self.x_];
                                
       end
       
       function [STOP_I] = stop_inner(self)
           STOP_I  = false;           
           %delta = norm(self.eval.cost_new - self.eval.cost_old_I);
           if self.iter_I > self.MAX_I_ITER
               STOP_I  = true;
           end            
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function plot(self)
           figure(1)
           cc=plot(self.x_(1:2:end),self.x_(2:2:end),'co')
           hold on
           axis([-0.1 0.4 -0.1 0.1 0 0.1])
           axis equal
           pause(0.1)
           delete(cc)
           axis([-0.1 0.4 -0.1 0.1 0 0.1])
           axis equal
       end
           
       
   end
end