function [ c_x ]= c_f(x,c,m,epsilon) 

Dfx = D_f(x,m,c);

if Dfx<0
    c_x = -Dfx+(1/2)*epsilon;    
elseif Dfx>=0 && Dfx<=epsilon
    c_x = (1/(2*epsilon))*(Dfx-epsilon)^2;
else
    c_x = 0;
end

end
