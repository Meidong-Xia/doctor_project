% The implement of new version of SVD required by HUAWEI
%
% Based on svd that comes with MATLAB
% 
% log:
%   - initialized by Meidong Xia on 01/15/2024
% 
%
% Rest of the code...
function [U,A,Vt]=cfp16_mysvd_4(mat)
    mat_half=hp_matrix_turn(mat);
    [U_single,A_single,Vt_single]=cfp32_svd_household(mat_half);
    U = hp_matrix_turn(double(U_single));
    A = hp_matrix_turn(double(A_single));
    Vt = hp_matrix_turn(double(Vt_single));
end