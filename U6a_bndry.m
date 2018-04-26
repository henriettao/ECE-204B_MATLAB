function [U] = U6a_bndry(t, n, ny)
U = -Inf*ones( n, n );
 
for j = 1:n
    for k = 1:n
        x = (j - 1)/(n - 1);
        y = (k - 1)/(n - 1);
        
        r = sqrt( (x - 0.5)^2 + (y - 0.5)^2 );
        
        if r < 0.1
           
            U(j, k) = 200-130*(heaviside(t-0.1));
            


        elseif r >= (n - 1)/(2*n);
            if y <= 0.5
                U(j, k) = 0;
            else
                U(j, k) = NaN;
            end
        end
    end
end
end