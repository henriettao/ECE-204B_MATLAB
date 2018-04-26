function [u] = u2b_bndry(t)
    u = [1, (t >=  2).*(3)];
end
