function [t, U_soln] = diffusion3d( kappa, h, U_init, U_bndry, t_int, nt )
	
    if ~isscalar(kappa) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument kappa is not a scalar' ) );
    end
    if ~isscalar(h) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument h is not a scalar' ) );
    end
    if ndims(U_init) ~= 3
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument U_init is not a three-dimensional matrix' ) );
    end
    if ~isa( U_bndry , 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument U_bndry is not a function handle' ) );
    end
    if ~all( size( t_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument t_int is not a 2-dimensional row vector' ) );
    end
    if ~isscalar( nt ) || ( nt ~= round( nt ) )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument nt is not a integer' ) );
    end
    
    delta_t = (t_int(2)-t_int(1))/(nt-1);
	r_check = (kappa*(delta_t)/(h))^2;
    
	if(r_check > (1/6))
    	n_t = ceil(((kappa*(t_int(2)-t_int(1)))/h)+1);
    	warning( '(c*delta_t/h)^2 returns %d >= 1, try nt = %d',r_check, n_t);
	end

    t0 = t_int(1);
    tf = t_int(2);
    dt = (tf - t0)/(nt - 1);
    t = linspace( t0, tf, nt );
 
    [nx, ny, nz] = size( U_init );
 
    U_soln = zeros( nx, ny, nz, nt );
    U_soln(:, :, :, 1) = U_init;
    
    r = kappa*dt/h^2;

    for it = 2:nt
        U_soln(:, :, :, it) = U_bndry( t(it), nx, ny, nz );
        
        for ix = 1:nx
            for iy = 1:ny
                for iz = 1:nz
                    if U_soln(ix, iy, iz, it) == -Inf
                        Utmp = U_soln(ix, iy, iz, it - 1);
                        U_soln(ix, iy, iz, it) = Utmp;

                        for dxyz = [-1, 1, 0, 0, 0, 0; 0, 0, -1, 1, 0, 0; 0, 0, 0, 0, -1, 1]
                            dix = ix + dxyz(1);
                            diy = iy + dxyz(2);
                            diz = iz + dxyz(3);

                            if ~isnan( U_soln(dix, diy, diz, it - 1) )
                                U_soln(ix, iy, iz, it) = U_soln(ix, iy, iz, it) + ...
                                    r*( U_soln(dix, diy, diz, it - 1) - Utmp );
                            end
                        end
                    end
                end
            end
        end
    end

end
