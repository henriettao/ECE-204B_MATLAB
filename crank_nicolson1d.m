% crank_nicolson1d
%
% the function can be used to solve insulated and non insulated boundary
% value problems

% We are solving the heat-conduction/diffusion equation on an interval
% a <= x <= b for a time interval t0 <= t <= tf. The initial conditions are
% given by the function u_init:[a, b] -> R and the boundary conditions are
% given by two functions a_bndry:[t0, t_f] -> R and
% b_bndry:[t0, tf] -> R.
%
% The solution U(x, t) must satisfy the initial and boundary conditions and
% must satisfy the heat-conduction/diffusion equation at every point in the
% region. We will approximate the solution on an nx x nt grid of points.
%
% Parameters
% ==========
% kappa The diffusion coefficient dictated by the system.
% x_rng The spacial interval [a, b] on which we will approximate the
% solution.
% t_rng The time interval [t0, tf] on which we will
% approximate the solution.
%
% u_init A real-valued function of a real variable describing the
% initial condition of the property at time t0nitial.
% u_bndry A vector-valued function returning a 2-dimensional column vector
% of the boundary values a_bndry(t) and b_bndry(t).
%
% nx We will approximate the solution at nx equally spaced points on [a, b].
% nt We will approximate the solution at nt equally spaced time steps for
% t0 <= t <= tf.
%
% Return Values
% =============
% x_out The vector x contains the nx points at which we are approximating the
% solution in the spatial dimension.
% t_out The vector t contains the nt points in time at which we will be
% approximating the solution.
% U_out The matrix U is an nx x nt matrix where U(k,ell) is the approximation of
% the solution the location x(k) and time t(ell). 

function [x_out, t_out, U_out] = crank_nicolson1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )


% Error Checking
	% ==============
	%
	% Checks if arguments are valid. Kappa, nx, and nt must be scalar.
	% x_rng and t_rng must be 2-dimensional row vectors. And u_init and
	% u_bndry must be functions
    % also ensures that nx is less than 4 
    % no need to check kappa*delta t/h*h crank_nicolson1d is implicit and
    % unconditionally stable
    % 
    h = (x_rng(2)- x_rng(1))/(nx-1);
	delta_t = (t_rng(2)-t_rng(1))/(nt-1);
	c = kappa*(delta_t)/(h*h);
    
	if ~isscalar(kappa)
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument kappa is not a scalar' ) );
	end
	if ~all( size( x_rng ) == [1, 2] )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument x_rng is not a 2-dimensional row vector' ) );
	end
	if ~isscalar( nx )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument nx is not a integer' ) );
	end
	if ~all( size( t_rng) == [1, 2] )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument t_rng is not a 2-dimensional row vector' ) );
	end
	if ~isscalar( nt ) || ( nt ~= round( nt ) )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument nt is not a integer' ) );
	end
	if ~isa( u_init, 'function_handle' )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument u_init is not a function handle' ) );
	end
	if ~isa( u_bndry , 'function_handle' )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument u_bndry is not a function handle' ) );
    end
    
    if(nx < 4)
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'nx is less than 4') );
    end
     if(c >= 0.5)
        warning( 'the arguments of %d and %d are sub-optimal',c,kappa);
     end
  

	% Initialization
	% ==============
	%
	%   Initializes the value of h, delta_t, c, x_out, t_out, and the
	%   matrix. In the special case where c >= 0.5 it is invalid. An
	%   appropriate nt is suggested to the user. The x_out and t_out are
	%   created by setting the evenly spaced intervals on x_rng and t_rng
	%   with nx and nt points. Matrix initializd by using the x_out and t_out
	%   vectors to create the first column and the first and last rows.
    %   creates the diagnoal matrix that is used to solve the system of
    %   linear equations 
    %   initializes the matrix which holds the solution to the diffusion
    %   equation
    %   initializes u_vector diff_vector and bndry_vector
    %   vectors and matrices are initalized to have the same dimersions 
    
	x_out= linspace(x_rng(1), x_rng(2), nx)';
	t_out= linspace(t_rng(1), t_rng(2), nt);
    
	matrix = zeros(nx,nt);
	matrix(:,1) = u_init(x_out(:,1));
    
	t_values = u_bndry(t_out);
	matrix(1,2:end) = t_values(1,1:end-1);
	matrix(nx,2:end) = t_values(2,1:end-1);

	u_vector = zeros(nx-2, 1);
	diff_vector = zeros(nx-2, 1);
	bndry_vector = zeros(nx-2, 1);
    
    
	r = diag((2*(1+c))*ones(nx-2, 1));
	r_top = diag(-c*ones((nx-3),1),-1);
	r_bottom = diag(-c*ones((nx-3),1),1);
	r_matrix = r_top+r+r_bottom;
    
    % Solving
    % =======
    %

    %   creates a system of linear equations to solve the values in matrix
    %   initalizes the column vector u_vector to the current column of the
    %   matrix 
    %   loops through the matrix and assignes the solution to the linear
    %   equation to the next column of the matrix
    %   U_out is then assigned as the matrix.
    
    %   Dirichlet boundary condition
    %   ==================
    %   if the boundry value problem is not insulated i.e a_bndry and b_ndry are numbers 
    %   the first and last row of matrix is assigned to the boundary
    %   considiotns 
    %   the missing values is calculated using the solution to the linear
    %   the matrix is assigned to U_out 
    %
    %   Insulated boundary 
    %   =================
    %   if the boundary value problem is insulated 
    %   that is one of the boundary conditions is set to Nan
    %   the matrix r_matrix ia altered to compensate for these changes 
    %   the missing values in the matrix is also assigned to the soltuion of
    %   the linear equation as before 
    %   and the boundary condition that was previous set to Nan is replace
    %   by using the surrounding values 
    % the matrix is assigned to U_out
   
	for j = 1:(nt-1)
   	for i = 2:(nx - 1)
       	u_vector(i-1) = matrix (i, j);
       	diff_vector(i-1) = matrix((i-1), j) - 2*matrix(i, j) + matrix((i+1), j);
   	end

   	if (isnan(t_values(2,j+1)) == 1)
        	r_matrix(1,1) = 2 + (2/3)*c;
        	r_matrix(1,2) = -(2/3)*c;
        	r_matrix(end,end) = 2 + (2/3)*c;
        	r_matrix(end,end-1) = -(2/3)*c;
        	known_vector = (2*u_vector) + (c*diff_vector);
        	matrix(:, j+1) = matrix(:, j+1)+[0; r_matrix\known_vector; 0];
        	matrix(end, j+1) = (4*matrix(end-2, j+1)) - (3*matrix(end-1, j+1));
   	else
        	bndry_vector(1) = c*t_values(1);
        	bndry_vector(end) = c*t_values(2);
        	known_vector = (2*u_vector) + (c*diff_vector) + bndry_vector;
        	matrix(:, j+1) = matrix(:, j+1)+[0; r_matrix\known_vector; 0];
   	end
  	 
	end
    
	U_out = matrix;

end