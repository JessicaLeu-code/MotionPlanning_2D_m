%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFS class file for xy-planning 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef CFS
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
       fun
       solver = 'CFS';
         
       % evaluation setting
       eval EVAL
       
       % calculation statistics
       iter_I = 1;
       iter_O = 1;
       total_iter = 0;
       
   end
   methods
       % _init_
       function self = CFS(val,val2,varargin)
                % get problem info
                self.obs = val;
                self.sys_info = val2;                
                self.nn = self.sys_info.H*self.sys_info.dim;                
                % initialize solution
                self.u = zeros(self.nn,1);                
                % initialize evaluation functions
                self.eval = EVAL(val2);
                if ~isempty(varargin)
                    self.x_ = varargin{1};
                else
                    self = self.get_path();
                end
                self.solver = self.sys_info.solver;
                self.fun = @(u) 0.5*u'*self.sys_info.Qaug*u + self.sys_info.paug'*u+self.sys_info.caug;
                                               
       end
       
       % main function
       function self = optimizer(self)
           self.iter_I = 1;
           self.iter_O = 1;
           self.u = zeros(self.nn,1);
           self.eval.cost_old = 10000;
           self.eval.cost_new = self.eval.get_cost(self.u);
           while self.eval.stop_outer(self.iter_O)~=true
               self.eval.cost_old = self.eval.cost_new;
               % get CFS               
               self = self.get_Gamma_AB();
               % do PSG
               self = self.iter_CFS();
               % next step
               self.iter_O = self.iter_O+1;
               self = self.get_path();
               %self.plot();
           end           
       end
              
       %%%%%%%%%%%%%%%% solver functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function self = iter_CFS(self)
           self.eval.u_old = self.u;
           if self.solver =='CFS'
               [self.u,fval,exitflag,output,lambda]= quadprog(self.sys_info.Qaug,self.sys_info.paug,self.Ainq,self.binq,[],[],[],[],[],self.sys_info.options); %quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options)
           elseif self.solver == 'SQP'
               [self.u,fval,exitflag,output,lambda] = fmincon(self.fun,self.eval.u_old,self.Ainq,self.binq,[],[],[],[],[],self.sys_info.options);
           end
           self.total_iter = self.total_iter + output.iterations;
           self.eval.cost_new = self.eval.get_cost(self.u);
           self.eval = self.eval.store_result(self.u);
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
       
       %%%%%%%%%%%%%%%% plotting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
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