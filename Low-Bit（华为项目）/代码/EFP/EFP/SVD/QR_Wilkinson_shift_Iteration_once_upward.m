function [tri,U]=QR_Wilkinson_shift_Iteration_once_upward(tri,U,i_min,i_max,shift)
n = i_max+1;
ct = i_hp_sub(i_min,shift);
% ct = i_min-shift;
x = i_hp_mul(tri(n,n),tri(n,n));
% x=tri(n,n)^2;
y = i_hp_mul(tri(n-1,n-1),tri(n-1,n));
% y=tri(n-1,n-1)*tri(n-1,n);
for k = n-1:-1:i_min
    [c,s,r] = csr(x,y);
    if k < n-1
        tri(k+1,k+2) = r;
    end
    A = i_hp_matrix_mul(complex([c,-s;s,c]),complex([tri(k+1,k+1),0;tri(k,k+1),tri(k,k)]));
    A = real(A);
%     A=[c,-s;s,c]*[tri(k+1,k+1),0;tri(k,k+1),tri(k,k)];
    x = A(1,1);
    y = A(1,2);
    tri(k,k+1) = A(2,1);
    tri(k,k) = A(2,2);
    [U(:,k+1),U(:,k)] = updatecsvv(c,s,U(:,k+1),U(:,k));
    [c,s,r] = csr(x,y);
    tri(k+1,k+1) = r;
    if k ~= ct+1
        A = i_hp_matrix_mul(complex([tri(k,k+1),tri(k,k);0,tri(k-1,k)]),complex([c,s;-s,c]));
        A = real(A);
%         A=[tri(k,k+1),tri(k,k);0,tri(k-1,k)]*[c,s;-s,c];
        x = A(1,1);
        tri(k,k) = A(1,2);
        y = A(2,1);
        tri(k-1,k) = A(2,2);
    else
        A = i_hp_matrix_mul(complex([tri(ct+1,ct+2),tri(ct+1,ct+1)]),complex([c,s;-s,c]));
        A = real(A);
%         A=[tri(ct+1,ct+2),tri(ct+1,ct+1)]*[c,s;-s,c];
        tri(ct+1,ct+2) = A(1);
        tri(ct+1,ct+1) = A(2);
    end
end
end