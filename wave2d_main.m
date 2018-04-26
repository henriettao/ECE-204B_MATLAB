nx = 301;
ny = 101;
U6c_init = zeros( nx, ny );
dU6c_init = zeros( nx, ny );

for ix = 1:nx
    for iy = 1:ny
        U6c_init(ix, iy) = exp( -0.1*((ix - 26)^2 + (iy - 51)^2) );
    end
end
[t6c, U6c] = wave2d( 1, 1, U6c_init, dU6c_init, @U6c_bndry, [0, 350], 710 );
frames6c = animate( U6c );