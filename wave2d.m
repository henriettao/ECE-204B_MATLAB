function [t, U_soln] = wave2d( c, h, U_init, dU_init, U_bndry, t_int, nt )
    
    if ~isscalar(c) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a scalar' ) );
    end
    if ~isscalar(h) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument h is not a scalar' ) );
    end
    if h<=0 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument h must be a positive number' ) );
    end
    
    if ndims(U_init) ~= 2
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument U_init is not a two-dimensional matrix' ) );
    end
    if ndims(dU_init) ~= 2
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument dU_init is not a two-dimensional matrix' ) );
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
	r_check = (c*(delta_t)/(h))^2;
    
	if(r_check > (1/4))
    	n_t = ceil(((c*(t_int(2)-t_int(1)))/h)+1);
    	warning( '(c*delta_t/h)^2 returns %d >= 1, try nt = %d',r_check, n_t);
	end

    t0 = t_int(1);
    tf = t_int(2);
    dt = (tf - t0)/(nt - 1);
    t = linspace( t0, tf, nt );
 
    [nx, ny] = size( U_init );
 
    U_soln = zeros( nx, ny, nt );
    U_soln(:, :, 1) = U_init;
    
    r = (c*dt/h)^2;

    U_soln(:, :, 2) = U_soln(:, :, 1) + dt*dU_init;

    for it = 3:nt
        U_soln(:, :, it) = U_bndry( t(it), nx, ny );
        
        for ix = 1:nx
            for iy = 1:ny
                if U_soln(ix, iy, it) == -Inf
                    Utmp = U_soln(ix, iy, it - 1);
                    U_soln(ix, iy, it) = 2*Utmp - U_soln(ix, iy, it - 2);
                    
                    for dxy = [-1 1 0 0; 0 0 -1 1]
                        dix = ix + dxy(1);
                        diy = iy + dxy(2);

                        if ~isnan( U_soln(dix, diy, it - 1) )
                            U_soln(ix, iy, it) = U_soln(ix, iy, it) + ...
                                r*( U_soln(dix, diy, it - 1) - Utmp );
                        end
                    end
                end
            end
        end
    end

end
