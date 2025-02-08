function y = f_genLaprnd(m,n,mu,sigma)
% f_genLaprnd generates m by n random matrix whose elements satisfy i.i.d.
% Laplacian distribution with mu mean and standard variance sigma.
%
% Based on:
%   http://en.wikipedia.org./wiki/Laplace_distribution
%
% Inputs:
%   - m: the row dimension of desired matrix
%   - n: the column dimension of desired matrix
%   - mu: mean
%   - sigma: standard variance
%
% Outputs:
%   - y: the desired m by n random matrix
%
% Default mu = 0, sigma = 1
%
% log:
%   - initialized by Meidong Xia on 9/27/2023
% 
%
% Rest of the code... 

if nargin == 2
    mu = 0;
    sigma = 1;
end

if nargin == 3
    sigma = 1;
end

u = rand(m,n)-0.5;
b = sigma/sqrt(2);
y = mu + sign(u) .* ( pi - b.* log( exp(pi./b) + (2-2.*exp(pi./b)) .* abs(u) ) );

end