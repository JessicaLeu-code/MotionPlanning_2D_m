function [Ainq, binq] = feasible_set(obs,x_,sys_info)
% Ainq*u <= binq

    % get # obstacles
    nobs = sys_info.num_obs;
    H = sys_info.H;
    Aaug = sys_info.Aaug;
    Baug = sys_info.Baug;
    dim = sys_info.dim;
    nstate = sys_info.nstate;
    in_max = sys_info.in_max;
    x0 = x_(1:nstate);
    dt = sys_info.dt;
    
    Ainq = [];
    binq = [];
    
    % environmental constraints
    for i = 1:nobs 
        if ~isempty(obs{i})
            switch obs{i}.shape
                case 'circle'            
                    for j = 1:H
                        x = x_(j*nstate+1:j*nstate+2);
                        % g(x*)+dg(x*)'(Ax0+Bu -x*)<=0
                        gx = g_f(obs{i}.A,obs{i}.c+obs{i}.v*dt*j,obs{i}.m,x);
                        dgx = dg_f(obs{i}.A,obs{i}.c,x);

                        Ainq = [Ainq; dgx'*Baug((j-1)*nstate+1:(j-1)*nstate+2,:)];
                        binq = [binq; dgx'*x-gx-dgx'*Aaug((j-1)*nstate+1:(j-1)*nstate+2,:)*x0];
                    end            
                case 'rectangle'
                    for j = 1:H
                        x = x_(j*nstate+1:j*nstate+2);
                        poly = obs{i}.poly+obs{i}.v*ones(1,4)*dt*j;
                        [dgx,gx,d] = d2poly(x',poly'); dgx=dgx';

                        Ainq = [Ainq; dgx'*Baug((j-1)*nstate+1:(j-1)*nstate+2,:)];
                        binq = [binq; gx-obs{i}.epsilon-dgx'*Aaug((j-1)*nstate+1:(j-1)*nstate+2,:)*x0];

                    end
            end
        end
        
    end
     
    % input constraints   (uncomment in inter_example) 
    Ainq = [Ainq; kron(eye(H),[eye(dim);-eye(dim)])];
    binq = [binq; kron(ones(H*2,1),in_max)];
    
            
        
end
