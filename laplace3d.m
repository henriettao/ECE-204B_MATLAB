function [U_out] = laplace3d( U )
	[n_x, n_y, n_z] = size( U );
	U_out = U;

	% Step 2
	u_to_w = zeros( n_x, n_y, n_z );
	w_to_u = zeros( 3, n_x * n_y * n_z );
	m = 0;

	for ix = 1:n_x
    	for iy = 1:n_y
        	for iz = 1:n_z
            	if U(ix, iy, iz) == -Inf
                	m = m + 1;
                	u_to_w(ix, iy, iz) = m;
                	w_to_u(:, m) = [ix, iy, iz]';
            	end
        	end
    	end
	end

	% Create the sparse system of linear equations
	M = spalloc( m, m, 7*m );
	b = zeros( m, 1 );

	for k = 1:m
    	c  = w_to_u(:,k);
   	 
    	p = c + [-1 0 0]';
    	val = U(p(1), p(2), p(3));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2), p(3))) = M(k, u_to_w(p(1), p(2), p(3))) + 1;
    	end
   	 
    	p = c + [1 0 0]';
    	val = U(p(1), p(2), p(3));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2), p(3))) = M(k, u_to_w(p(1), p(2), p(3))) + 1;
    	end
   	 
    	p = c + [0 -1 0]';
    	val = U(p(1), p(2), p(3));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2), p(3))) = M(k, u_to_w(p(1), p(2), p(3))) + 1;
    	end
   	 
    	p = c + [0 1 0]';
    	val = U(p(1), p(2), p(3));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2), p(3))) = M(k, u_to_w(p(1), p(2), p(3))) + 1;
    	end
   	 
    	p = c + [0 0 -1]';
    	val = U(p(1), p(2), p(3));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2), p(3))) = M(k, u_to_w(p(1), p(2), p(3))) + 1;
    	end
   	 
    	p = c + [0 0 1]';
    	val = U(p(1), p(2), p(3));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2), p(3))) = M(k, u_to_w(p(1), p(2), p(3))) + 1;
    	end
	end

	w = M \ b;

	for k = 1:m
    	p  = w_to_u(:,k);
    	U_out(p(1), p(2), p(3)) = w(k);
	end
end