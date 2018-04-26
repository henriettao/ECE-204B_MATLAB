% REPLACE THIS TEMPLATE WITH YOUR COMMENTS AND SIGNATURE
% laplace2d
% Copy the description from 5.1 here
%
% Parameters
% ==========
%    U:		The input matrix of the problem over which initial and boundary conditions are placed.
%
% Return Values
% =============
%    U_out:	The approximation of the laplacean equation solution to the initial matrix.
% Error and Warning Checking
% ==========================
%
%   An exception should be thrown if the input matrix U, is not a two-dimensional matrix.
%  

%  Your commands here...

% Initialization
% ==============
%
%   Your description here.
% Stores the dimensions of U for future uses and assigns U as U_out.
% Initializes a  u_to_w matrix to a zero matrix of the same size as U.
% Initializes a w_to_u matrix to a zero matrix with two rows and length equal to the % % % number of -Inf in U.
% Initializes constant m to 0 in preparation of mapping points.
% Matrix M and column vector b is initialized in order to solve the system of linear equations.

%   Your commands here...

% Mapping the unknown points to a unique number from 1 to m
% =========================================================
%
%   Your description here.
% In the u_to_w matrix, an index number of each -Inf of U is stored in the same position % in u_to_w matrix as the -Inf in U.
% In the w_to_u matrix, the coordinate of each -Inf in U is stored as a column vector.

%     Your commands here...

% Creating and solving a system of linear equations
% =================================================
%
%   Your description here.
% Matrix M and column vector b is initialized in order to solve the system of linear equations.
% Approximate the solutions by iterating through all points and %checking the neighbouring points: to the left, to the right, %above, and below the point. If the value is an insulated %boundary (NaN), nothing is done. If the value is an Dirichlet %boundary, then 1 is subtracted from the kth index on the % %diagonal of M, and the value is subtracted from the kth index of b. If the value is an unknown point (-Inf), then 1 is subtracted from the kth index on the diagonal of M, and 1 is added to the coordinate position in M.
% Then the solutions is solved by getting w = M \ b.

% Your commands here...

% Substituting the values back into the matrix U_out
% ===================================================
%
%   Your description here.

% Your commands here...

function [U_out] = laplace2d( U )
	[n_x, n_y] = size( U );
	U_out = U;

	% Step 2
	u_to_w = zeros( n_x, n_y );
	w_to_u = zeros( 2, n_x * n_y );
	m = 0;

	for ix = 1:n_x
    	for iy = 1:n_y
        	if U(ix, iy) == -Inf
            	m = m + 1;
            	u_to_w(ix, iy) = m;
            	w_to_u(:, m) = [ix, iy]';
        	end
    	end
	end
	% Create the sparse system of linear equations
	M = spalloc( m, m, 5*m );
	b = zeros( m, 1 );

	for k = 1:m
    	c  = w_to_u(:,k);
   	 
    	p = c + [-1 0]';
    	val = U(p(1), p(2));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2))) = M(k, u_to_w(p(1), p(2))) + 1;
    	end
   	 
    	p = c + [1 0]';
    	val = U(p(1), p(2));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2))) = M(k, u_to_w(p(1), p(2))) + 1;
    	end
   	 
    	p = c + [0 -1]';
    	val = U(p(1), p(2));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2))) = M(k, u_to_w(p(1), p(2))) + 1;
    	end
   	 
    	p = c + [0 1]';
    	val = U(p(1), p(2));
    	if ~isnan(val) && (val ~= -Inf)
        	M(k, k) = M(k, k) - 1;
        	b(k) = b(k) - val;
    	elseif val == -Inf
        	M(k, k) = M(k, k) - 1;
        	M(k, u_to_w(p(1), p(2))) = M(k, u_to_w(p(1), p(2))) + 1;
    	end
	end
    
	w = M \ b;

	for k = 1:m
    	p  = w_to_u(:,k);
    	U_out(p(1), p(2)) = w(k);
    	% Copy the value from w into U_out
	end
end