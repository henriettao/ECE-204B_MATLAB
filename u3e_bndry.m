function [u] = u3e_bndry( t )
    nt = length( t );
    
    u = zeros( 2, nt );
    for k = 1:nt
        u(:, k) = [100; NaN];
    end
end