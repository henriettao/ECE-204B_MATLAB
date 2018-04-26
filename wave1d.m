% wave1d
%
% The function is used to approximate wave equations for given parameters.
%
% Parameters
% ==========
% c: 		The diffusion coefficient dictated by the system.
% x_rng: 	The two dimensional row vector that describes the interval of the
%  			boundary in space
% t_rng: 	The two dimensional row vector that describes the interval of the
%  			boundary in time
%
% u_init: 	The function where the x-value is evaluated at a_init and b_init
% u_bndry: 	The function where a_bndry and b_bndry is evaluated at all t values
% 			and returns a two dimensional column vector
%
% nx: 		The number of points equally spaced on the interval [a,b]
% nt: 		The number of points equally spaced on the interval [t0, tfinal]
%
% Return Values
% =============
% x_out: 	The vector x is a vector of the nx points across interval [a, b]
% t_out:	The vector t is a vector of the nt points across interval 
%			[t0, tfinal]
% U_out:	The matrix U is a matrix of values approximating the solution of
%			the wave equation of size nx*nt
 
function [x_out, t_out, U_out] = wave1d( c, x_rng,  nx, t_rng, nt, u_init, du_init, u_bndry )
 
	% Error Checking
	% ==============
	%
	% Checks if arguments are valid. c must be scalar. nt and nx must
	% be integers. x_rng and t_rng must be 2-dimensional row vectors. u_init, 
	% du_init, and u_bndry must be functions.
	% Shows warning if (c*delta_t/h)^2 is greater than 1 and displays an 
	% appropriate nt to use instead.
	%
    
	if ~isscalar(c)
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument kappa is not a scalar' ) );
	end
	if ~all( size( x_rng ) == [1, 2] )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument x_rng is not a 2-dimensional row vector' ) );
	end
	if ~isscalar( nx ) || ( nx ~= round( nx ) )
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
	if ~isa( du_init, 'function_handle' )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument u_init is not a function handle' ) );
	end
	if ~isa( u_bndry , 'function_handle' )
    	throw( MException( 'MATLAB:invalid_argument', ...
    	'the argument u_bndry is not a function handle' ) );
	end

	h_check = (x_rng(2)- x_rng(1))/(nx-1);
	delta_t_check = (t_rng(2)-t_rng(1))/(nt-1);
	r_check = (c*(delta_t_check)/(h_check))^2;
    
	if(r_check > 1)
    	n_t = ceil(((c*(t_rng(2)-t_rng(1)))/h_check)+1);
    	warning( '(c*delta_t/h)^2 returns %d >= 1, try nt = %d',r_check, n_t);
	end
 
	% Initialization
	% ==============
	%
	%   Initializes the value of h, delta_t, r, x_out, t_out, and the matrix.
	%	The x_out and t_out are created by setting the evenly spaced intervals
	%	on x_rng and t_rng with nx and nt points. Matrix initialized by using 
	%	the x_out and t_out vectors to create the first column and the first 
	%	and last rows using the u_init, du_init, and u_bndry functions.
    
	h = (x_rng(2)- x_rng(1))/(nx-1);
	delta_t = (t_rng(2)-t_rng(1))/(nt-1);
	r = (c*(delta_t)/(h))^2;
    
	x_out= linspace(x_rng(1), x_rng(2), nx)';
	t_out= linspace(t_rng(1), t_rng(2), nt);
    
	matrix = zeros(nx,nt);
	matrix(:,1) = u_init(x_out(:,1));
    
	t_values = u_bndry(t_out);
	matrix(1,2:end) = t_values(1,1:end-1);
	matrix(nx,2:end) = t_values(2,1:end-1);
    
	% Solving
	% =======
	%
	%   Calculates values for the rest of the matrix. Solves first column of 
	%	matrix using euler's method. Solves the rest of the matrix using 
	%	expression. The matrix is assigned to U_out.
	%
	%   Insulated boundary
	%   =================
	%   If the boundary value problem is insulated, that means that one of the
	%	boundary conditions is set to NaN. Every column is checked to see that
	%	the a or b boundary is not NaN. If it is NaN then an appropriate number
	%	is calculated using surrounding values and replaces the boundary.
	%   
   
	for j = 1:(nt-1)
    	for i = 2:(nx - 1)
        	if (j == 1)
            	matrix(i,j+1) = matrix(i, j) + (delta_t*du_init(i)) + ((r/2)*(matrix(i-1,j) - (2*matrix(i,j)) + matrix(i+1,j)));
        	else
            	matrix(i,j+1) = (2*matrix(i,j)) - matrix(i,j-1) + (r*(matrix(i-1,j) - (2*matrix(i,j)) + matrix(i+1,j)));
        	end
    	end
 
    	if (isnan(t_values(1,j+1)) == 1)
        	matrix(1, j+1) = ((4/3)*matrix(2, j+1)) - ((1/3)*matrix(3, j+1));
    	end
    	if (isnan(t_values(2,j+1)) == 1)
        	matrix(end, j+1) = ((4/3)*matrix(end-1, j+1)) - ((1/3)*matrix(end-2, j+1));
    	end
	 
	end
    
	U_out = matrix;
 
end