function [ dc_x ]= dc_f(x,c,m,epsilon) 

Dfx = D_f(x,m,c);
dDfx = dD_f(x,c);

if Dfx<0
    dc_x = -dDfx;    
elseif Dfx>=0 && Dfx<=epsilon
    dc_x = (1/epsilon)*(Dfx-epsilon)*dDfx;
else
    dc_x = [0;0];
end

end
