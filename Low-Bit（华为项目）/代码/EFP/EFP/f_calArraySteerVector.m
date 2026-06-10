function a = f_calArraySteerVector(M,theta)
% f_calArraySteerVector generates array steer vector given the angle and size.
%
% Inputs:
%   - M: the size of array steer vector
%   - theta: the angle of array steer vector
%
% Outputs:
%   - a: the array steer column vector
%
% log:
%   - initialized by Meidong Xia on 9/26/2023
% 
%
% Rest of the code... 

% idx = (1:M)-M/2;
idx = 0:M-1;
a = exp(1i*pi.*idx*sin(theta)).'/sqrt(M);
end