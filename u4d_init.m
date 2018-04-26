function [u] = u4d_init(x)
% 	u = sin(x);
    u = (1*x.*(heaviside(x) - heaviside(x-1))) + ((-1/3)*x + (4/3)).*(heaviside(x-1) - heaviside(x-4));
    
    
end