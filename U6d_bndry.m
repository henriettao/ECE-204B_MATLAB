function [U] = U6d_bndry( t, nx, ny, nz )
    U = -Inf*ones( nx, ny, nz );
    U(end,:,:) = sin(t).*(t <= 2*pi);

    
    U(:,[1,end],:) = NaN;
    U(:,:,[1,end]) = NaN;
    
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                x = (ix - 1)/(nx - 1);
                y = (iy - 1)/(ny - 1);
                z = (iz - 1)/(nz - 1);
                
                % Modify this to determine which points at the end
                % of the region will be set to NaN to reflect and
                % focus the signal.  Consider the example where 
                % the points were determined by the distance from 
                % the centre of a circle.
        
%                 if x == 0
%                     U(ix, iy, iz) = NaN;
%                 end
   `            if sqrt((y-0.25)^2+ (z-0.5)^2)== x
                 U(ix,iy,iz)=NaN;
                end
                
            end
        end
    end
end
