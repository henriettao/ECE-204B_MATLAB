function [U] = U6c_bndry( t, nx, ny )
    U = -Inf*ones( nx, ny );
 
    for ix = 1:nx
        for iy = 1:ny
            % (x, y) is a point in [0, 1] x [0, 1]
        
            % Determine if a point is outside a circle of radius 0.5
            % centred at the point (0.5, 0.5)
            if sqrt( (ix - 101)^2 + (iy - 51)^2 ) >= 100
                U(ix, iy) = NaN;
            end
        end
    end
end