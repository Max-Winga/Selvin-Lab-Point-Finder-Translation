function rectang_func = rectang(p,rr)
% p : how many fringes

tt = floor(length(rr)/(2*p));
r_pi = pi*ones(tt,tt);
r_0 = zeros(tt,tt);
fringe = [r_pi,r_0;r_0,r_pi];
rectang_func = repmat(fringe,p);
end