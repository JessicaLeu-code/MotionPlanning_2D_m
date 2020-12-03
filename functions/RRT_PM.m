%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RRT class file for mobile robots (xy planning) 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef RRT_PM
   properties
       % setup parameters
       obs cell
       sys_info struct
       goal {mustBeNumeric}
       region_g {mustBeNumeric}                    % goal region size
       region_s {mustBeNumeric}                    % samplign region size
       sample_off {mustBeNumeric}                  % sampling region center 
       SOLVER = 'RRT'                             % 'RTT*' or 'RRT'
       
       % intermediate states
       parent {mustBeNumeric}
       newNode {mustBeNumeric}
       isGoalReached logical
       toNode_dis {mustBeNumeric}                  % vector
       
       % outputs
       all_nodes {mustBeNumeric}                   % all nodes
       total_dis {mustBeNumeric}                   % distance away from initial point
       route {mustBeNumeric}                       % output route
       fail = 0;
           
       % sover settings
       MAX_ITER = 100;                            % maximum number of samples
       bi = 0.7;                                   % probability of 'not' sample the goal       
       rrt_stepSize = 0.015;                     
       rrt_rearrange_range = 0.03;       
             
       % calculation statistics
       node_num = 1;              
       
   end
      
   methods
       %%% _init_
       function self = RRT_PM(val,val2,val3,val4,val5,varargin)           
                % get problem info
                self.obs = val;
                self.sys_info = val2;                
                self.goal = self.sys_info.goal_th;
                self.region_g = val3;
                self.region_s = val4;
                self.sample_off = val5;
                if ~isempty(varargin)                    
                    self.SOLVER = varargin{1};
                end   
       end
              
       %%% Main function
       function self = RRTfind(self)
           % initialize tree
           self.newNode = self.sys_info.x0;
           self.all_nodes = [-1; self.newNode];
           self.total_dis = [0];
           self = self.goal_reached();
           % generate rrt
           switch self.SOLVER
               case 'RRT'
                   while self.isGoalReached ~=true
                       self = self.getNode();
                       self = self.addNode();               
                       self = self.goal_reached();               
                   end
               case 'RRT*'
                   while self.isGoalReached ~=true
                       self = self.getNode();
                       self = self.addNode();
                       self = self.arrangeNode();               
                       self = self.goal_reached();               
                   end
           end
           % find route
           self.route = self.newNode;
           while self.parent ~= -1
               self.route = [ self.all_nodes(2:end,self.parent) self.route];
               self.parent = self.all_nodes(1,self.parent);
           end
       end
       %%%%%%%%%%%%%%%% supporting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% get feasible node
       function  self = getNode(self)
           self = self.getRandNode();
           [isFeasible, self]= self.feasible();
           while isFeasible ~=true
               self = self.getRandNode();
               [isFeasible, self]= self.feasible();               
           end    
           
       end               
       %%%% get random node
       function self = getRandNode(self)
           pp = rand;
           if pp<self.bi
               % specified sample range 
               sampleNode = (rand(self.sys_info.nstate,1)-0.5).*self.region_s*2+self.sample_off;
           else
               sampleNode = self.sys_info.goal_th;
           end

           self.parent = 1;
           dis = norm(self.all_nodes(2:end,1)-sampleNode);           
           for i = 2:self.node_num
               
               dis_ = norm(self.all_nodes(2:end,i)-sampleNode);
               if dis_<dis
                   dis = dis_;
                   self.parent = i;
               end
           end
           self.toNode_dis(self.parent) = self.rrt_stepSize;
           self.newNode = self.all_nodes(2:end,self.parent)+(sampleNode-self.all_nodes(2:end,self.parent))*self.rrt_stepSize/norm(self.all_nodes(2:end,self.parent)-sampleNode);
        end
       
       %%% rearrange parent-children pares in RRT*
       function self = arrangeNode(self)
           % get distance to new node
           self.toNode_dis = [];
           for i = 1:self.node_num -1
               dis_ = norm(self.all_nodes(2:end,i)-self.newNode);               
               self.toNode_dis = [self.toNode_dis dis_];               
           end
           % fund nearby nodes
           ids = find(self.toNode_dis <self.rrt_rearrange_range);
           % initialize rearrangement
           parent_candidate = self.parent;
           dis_ = self.total_dis(end);
           % find new arrangement
           for i = 1:length(ids)
               % find new parents for previous nodes
               if self.total_dis(ids(i))>(self.total_dis(end)+self.toNode_dis(ids(i)))
                   self.all_nodes(1,ids(i)) = self.node_num; 
                   self.total_dis(ids(i)) = (self.total_dis(end)+self.toNode_dis(ids(i)));
               end
               % find new parents for the newest nodes
               if dis_ >(self.total_dis(ids(i))+self.toNode_dis(ids(i)))
                   parent_candidate = ids(i);
                   dis_ = (self.total_dis(ids(i))+self.toNode_dis(ids(i)));
               end      
           end
           % save change for newest node
           self.parent = parent_candidate;
           self.all_nodes(1,end) = self.parent;
           self.total_dis(end) = dis_;
       end
       
       %%%%%%%%%%%%%%%%%%% supplementary functions %%%%%%%%%%%%%%%%%%%%%%%%     
       function [isFeasible, self]= feasible(self)
           for j = 1:self.sys_info.num_obs
               if ~isempty(self.obs{j})
                   switch self.obs{j}.shape
                       case 'circle' 
                           if norm(self.newNode - self.obs{j}.c)<self.obs{j}.m
                              isFeasible = false;
                              break;
                           else
                              isFeasible = true;
                           end          
                       case 'rectangle'
                           region_obs = [self.obs{j}.w/2;self.obs{j}.l/2];
                           if (self.obs{j}.c-region_obs)<self.newNode
                               if self.newNode<(self.obs{j}.c+region_obs)
                                   isFeasible = false;
                                   break;
                               else
                                   isFeasible = true;
                               end
                           else
                               isFeasible = true;
                           end 
                    end
                   
               end
           end
       end
       
       function self = addNode(self)           
           self.all_nodes = [self.all_nodes [self.parent; self.newNode]];
           self.total_dis = [self.total_dis (self.total_dis(self.parent) + self.toNode_dis(self.parent))];
           self.node_num = self.node_num + 1;
           
       end
             
       function self = goal_reached(self)
           %self = self.get_ee();
           self.isGoalReached = false;           
           if (self.goal-self.region_g)<self.newNode
               if self.newNode<(self.goal+self.region_g)
                   self.isGoalReached = true;
                   self.newNode = self.sys_info.goal_th;
                   self.parent = self.node_num;
                   self.toNode_dis(self.parent) = norm(self.all_nodes(2:end,self.parent)-self.newNode);
                   self = self.addNode();
               end
           end
           
           if self.node_num>self. MAX_ITER
               disp('Failed to find path.')
               self.fail = true;
               self.isGoalReached = true;
           end
                      
       end
       
     %%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      function self = plotNode(self)
%            %figure(self.sys_info.fighandle)
%            %plot3(self.newNodeXYZ(1),self.newNodeXYZ(2),self.newNodeXYZ(3),'*r')
%            %pause(0.01)
%      end
           
       
   end
end