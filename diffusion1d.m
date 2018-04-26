%   diffusion1d
%
%   The equation depends on multiple variables, so the partial derivatives 
%   need to be dealt with as such. This means that the equation needs to be
%   solved per variable to get the final answer.
%
% Parameters
% ==========
%    kappa
%      kappa: is the thermal diffusivity of the material 
%    x_rng
%      x_rng: is a two dimensional row vector that describes the interval of the 
%      boundary in space
%    t_rng
%      t_rng: is a two dimensional row vector that descibes the interval of the 
%      boundary condition in time  
%    u_init
%       u_int: uint(x) is a function where the x-value is evaluated at ainit and binit
%    u_bndry
%      u_bndry: ubndry(t) is a function where the value at abndry and bbndry is evaluated at t
%    nx
%      nx: the number of points equally spaced on the interval [a,b]
%    nt
%      nt: the number of points equally spaced on the interval [t0, tfinal]
% Return Values
% =============
%    x_out
%      x_out: The vector x is a vector of the nx points across interval [a, b]
%    t_out
%      t_out: The vector t is a vector of the nt points across interval [t0, tfinal]
%    U_out
%      U_out:The matrix U is a vector of values approximating the solution of size nx*nt
function [x_out, t_out, U_out] = diffusion1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )
    % Error Checking
    % ==============
    %
    % Checking if arguments are valid. Kappa, nx, and nt must be scalar. 
    % x_rng and t_rng must be 2-dimensional row vectors. And u_init and 
    % u_bndry must be functions. 
    
    if ~isscalar(kappa) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument kappa is not a scalar' ) );
    end
    if ~all( size( x_rng ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_rng is not a 2-dimensional row vector' ) );
    end
    if ~isscalar( nx) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument nx is not a scalar' ) );
    end
    if ~all( size( t_rng) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument t_rng is not a 2-dimensional row vector' ) );
    end
    if ~isscalar( nt) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument nt is not a scalar' ) );
    end
    if ~isa( u_init, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_init is not a function handle' ) );
    end
    if ~isa( u_bndry , 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a function handle' ) );
    end
    
    % Initialization
    % ==============
    %
    %	Initializes the value of h, delta_t, c, x_out, t_out, and the 
    %   matrix. In the special case where c >= 0.5 it is invalid. An 
    %   appropriate nt is suggested to the user. The x_out and t_out are 
    %   created by setting the evenly spaced intervals on x_rng and t_rng 
    %   with nx and nt points. Matrix initializd by using the x_out and t_out
    %	vectors to create the first column and the first and last rows.
    
	h = (x_rng(2)- x_rng(1))/(nx-1);
    delta_t = (t_rng(2)-t_rng(1))/(nt-1);
    c = kappa*(delta_t)/(h*h);

    if(c >= 0.5)
        n_t = ceil(((kappa* (t_rng(2)-t_rng(1))) / (0.5*(h*h)))+1);
        throw( MException( 'MATLAB:invalid_argument', ...
        'kappa*dt/h^2 = %d is >= 0.5, try using %d instead for nt', c, n_t) );
    end
    
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
    %   Loops through the matrix and uses the finite-difference equation to
    %   comput and assign the clues to the appropriate spot in the matrix.
    %   U_out is then assigned as the matrix.

	for j = 1:(nt-1)
       for i = 2:(nx - 1)
           matrix(i, (j + 1)) = (matrix(i, j)) + c*((matrix((i- 1), j)) - (2*(matrix(i, j))) + (matrix((i + 1), j)));
       end
    end

   U_out = matrix;
    
end
