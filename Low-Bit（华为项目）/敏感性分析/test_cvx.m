tic
idx = 0;
f = 0;
t = parpool(12);
parfor kk = 1:100
    k  = 3;
    P = 3;
    p = f_function(k,P);
    
    idx = idx+1;
    f = f+sum(p);
end
delete(t);
f./idx
toc

function p = f_function(k,P)
cvx_begin quiet
    variable p(k) nonnegative
    variable t
    minimize(t)
    subject to
        for i = 1:k
            inv_pos(p(i)) <= t
        end
        sum(p)==P
        p(1)>=2
cvx_end
end