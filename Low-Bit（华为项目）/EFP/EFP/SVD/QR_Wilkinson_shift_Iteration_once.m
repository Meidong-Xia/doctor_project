function [tri,U,V]=QR_Wilkinson_shift_Iteration_once(tri,U,V,i_min,i_max,shift)
n = i_max+1;
tri2 = real(i_hp_mul(tri(i_min+1,i_min+1),tri(i_min+1,i_min+1)));
x = real(i_hp_sub(tri2,shift));
% x = tri(i_min+1,i_min+1)*tri(i_min+1,i_min+1)-shift;
y = real(i_hp_mul(tri(i_min+1,i_min+1),tri(i_min+1,i_min+2)));
% y=tri(i_min+1,i_min+1)*tri(i_min+1,i_min+2);
for k = i_min:n-2
    [c,s,r] = csr(x,y);
    A = real(i_hp_matrix_mul(complex([tri(k+1,k+1),tri(k+1,k+2);0,tri(k+2,k+2)]),complex([c,s;-s,c])));
%     A=[tri(k+1,k+1),tri(k+1,k+2);0,tri(k+2,k+2)]*[c,s;-s,c];
    x = A(1,1);
    tri(k+1,k+2) = A(1,2);
    y = A(2,1);
    tri(k+2,k+2) = A(2,2);
    [V(:,k+1),V(:,k+2)] = updatecsvv(c,s,V(:,k+1),V(:,k+2));% 低精度
    if(k>i_min)
        tri(k,k+1) = r;
    end
    [c,s,r] = csr(x,y);
    tri(k+1,k+1) = r;
    [U(:,k+1),U(:,k+2)] = updatecsvv(c,s,U(:,k+1),U(:,k+2));% 高精度
    if k ~= n-2
        A = real(i_hp_matrix_mul(complex([c,-s;s,c]),complex([tri(k+1,k+2),0;tri(k+2,k+2),tri(k+2,k+3)])));
%         A = [c,-s;s,c]*[tri(k+1,k+2),0;tri(k+2,k+2),tri(k+2,k+3)];
        x = A(1,1);
        y = A(1,2);
        tri(k+2,k+2) = A(2,1);
        tri(k+2,k+3) = A(2,2);
    else
        A = real(i_hp_matrix_mul(complex([c,-s;s,c]),complex([tri(n-1,n);tri(n,n)])));
%         A=[c,-s;s,c]*[tri(n-1,n);tri(n,n)];
        tri(n-1,n) = A(1);
        tri(n,n) = A(2);
    end
end
end