% Physical parameters
diff_coeff = 0.01;
c_bath =  1.1198;
T_bath  = 6.03;
T_final = 12;
 
% Numerical parameters
n = 13;
nt = 251;
 
U_init = zeros( n, n, n );
U_init(:, :, [1, end]) = c_bath;
 
for i = 1:n
    for j = 1:n
        for k = 1:n
            x = (i - 1)/(n - 1);
            y = (j - 1)/(n - 1);
        
            r = sqrt( (x - 0.5)^2 + (y - 0.5)^2 );
        
            if r >= 0.5
                U_init(i, j, k) = c_bath;
            end
        end
    end
end

% Define U6b_tmp to call U6b_bndry with the two last parameters
U6b_tmp = @(t, n1, n2, n3)(U6b_bndry( t, n1, n2, n3, c_bath, T_bath ));

% Solve for the diffusion
[t, U_out] = diffusion3d( diff_coeff, 1/(n - 1), U_init, U6b_tmp, ...
                          [0, T_final], nt );
 
% Find the minimum and maximum concentrations for each point in time
c_max = zeros( 1, nt );
c_min = zeros( 1, nt );
 
for k = 1:nt
    c_max(k) = max( max( max( U_out(:, :, :, k) ) ) );
    c_min(k) = min( min( min( U_out(:, :, :, k) ) ) );
end
 
% Plot the minimum and maximum temperatures
plot( t, c_max );
hold on;
plot( t, c_min );
[c_min(end), c_max(end)]