function [x_out, u_out] = bvp( c, g, x_bndry, u_bndry, n )
% bvp
% approximates the general solution to a LODE with boundary values 
%
% Parameters
% ==========
%    c 
%  % vector of real constant values
%    g
%  % the function g(x) is the forced response of the system 
%    x_bndry
%    u_bndry
% % x_bndry and u_bndry describes the boundaries of the boundary value problem 
%    n
%  % The interval is divided into into n intervals to get a better approximation.
% Return Values
% =============
%    x_out
%  % The vector x_out is a column vector corresponding to the points between the upper bound and lower bound
%    u_out
%  % The vector u_out is a column vector corresponding to the function evaluated at the x boundaries

% Argument Checking
% Checks to see of c is a 3-dimensional column vector
    if ~all( size( c ) == [3, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a 3-dimensional column vector' ) );
    end
% Checks to see of g is a function
    if ~isa( g, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument g is not a function handle' ) );
    end
% Checks to see of x_bndry is a 3-dimensional column vector
    if ~all( size( x_bndry ) == [2, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_bndry is not a 2-dimensional column vector' ) );
    end
% Checks to see of u_bndry is a 3-dimensional column vector
    if ~all( size( u_bndry ) == [2, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a 2-dimensional column vector' ) );
    end
% Checks to see of n is a scalar
    if ~isscalar( n ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n is not a scalar' ) );
    end

% Sets x_out as all the points at n intervals between the x_boundry
x_out = linspace(x_bndry(1), x_bndry(2), n);

% Sets h as the distance bewteen each point on the intervals
h = (x_bndry(2)-x_bndry(1))/(n-1);

% Creates a super diagonal with d-, d, and d+ to put in a matrix
d = diag((2*h*h*c(3)-4*c(1))*ones((n-2),1));
d_minus = diag((2*c(1) - h*c(2))*ones((n-3),1),-1);
d_plus = diag((2*c(1)+h*c(2))*ones((n-3),1),1);
d_matrix = d_minus+d+d_plus;

% Sets x_vals as the intervals between the x_bndry
x_vals = linspace(x_out(2), x_out(end-1), n-2);

% Calculates b and sets toe b vector at each point on the interval
b = 2*h*h*g(x_vals);
b(1) = b(1) - ((2*c(1) - h*c(2))*u_bndry(1));
b(end) = b(end) - ((2*c(1)+h*c(2))*u_bndry(2));

% Solves the system of linear equations and concatentates the n points to the result
u_out = [u_bndry(1) (d_matrix\b')' u_bndry(2)];
end