nx = 25;
ny = 25;
nz = 25;
nt = 151;
U6d_init = zeros( nx, ny, nz );
dU6d_init = zeros( nx, ny, nz );
[t6d, U6d] = wave3d( 0.05, 1/(nx - 1), U6d_init, dU6d_init, ...
                     @U6d_bndry, [0, 40], nt );
frames6d = animate( U6d );