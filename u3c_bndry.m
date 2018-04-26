function [u] = u3c_bndry(t)
    nt = length( t );
    
    u = zeros( 2, nt );
    
    for k = 1:nt
        if t(k) < 5
            u(:, k) = [0; 2];
        else
            u(:, k) = [0; NaN];
        end
    end
end