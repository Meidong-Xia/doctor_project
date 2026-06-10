function [U,A,V,err]=cfp16_svd_lanczos(mat)
[size_matm,size_matn] = size(mat);
if size_matm > size_matn
    islong = 1;
else
    islong = 0;
end
if islong
    mat = mat';
end
[U,A,V] = newbilan(mat);
T = U;
A = A';
U = V;
V = T;
err = norm(U*A*V'-mat',"fro")/norm(mat,"fro");
[U,A,Vt] = bi_diag_svd(U,A,V);
[m0,n0] = size(A);
[lenu,~] = size(U);
[lenv,~] = size(Vt);
n = min(m0,n0);
for i = 1:n-1
    if real(A(i,i)) < 0
        A(i,i) = -A(i,i);
        for k = 1:lenu
            U(k,i) = -U(k,i);
        end
    end
end
for i = 1:n-1
    for j = i:n
        if(abs(A(i,i)) < abs(A(j,j)))
            tmp1 = A(i,i);
            A(i,i) = A(j,j);
            A(j,j) = tmp1;
            for k = 1:lenu
                temp2 = U(k,i);
                U(k,i) = U(k,j);
                U(k,j) = temp2;
            end
            for k = 1:lenv
                tmp2 = Vt(i,k);
                Vt(i,k) = Vt(j,k);
                Vt(j,k) = tmp2;
            end
        end
    end
end
if ~islong
    T = Vt';
    Vt = U';
    U = T;
    A = A';
end
V = Vt';

for i=1:min(lenu,lenv)
    if(real(U(1,i))<0)
        U(:,i)=-U(:,i);
        V(:,i)=-V(:,i);
    end

end

end























