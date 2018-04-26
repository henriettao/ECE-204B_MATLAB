function [u] = u4b_bndry(t)
    u = [(1 - cos(t)).*(t <= 2*pi); NaN*t];
end
