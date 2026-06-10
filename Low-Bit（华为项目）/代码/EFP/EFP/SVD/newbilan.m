function [U,B,V] = newbilan(mat)
% mat=mat';
% mat = H;
[m,n] = size(mat);
if m > n
    mat = mat';
    isT = 1;
else
    isT = 0;
end
[m,n] = size(mat);
mini = min(m,n);
p0 = zeros(m,1);
p0(1) = 1;
beta = zeros(mini+1);
alpha = zeros(mini);
U = zeros(m,mini+1);
V = zeros(n,mini);
beta(1) = 1;
U(:,1) = p0;
B = zeros(mini+1,mini);
r = zeros(n,mini+1);
P = zeros(m,mini);
for j = 1:mini
    if j ~= 1
        mu = i_hp_matrix_mul(complex(mat'),complex(U(:,j)));
        bv = i_hp_matrix_mul_e(complex(V(:,j-1)),complex(beta(j)));
        r(:,j) = i_hp_matrix_sub(complex(mu),complex(bv));
%         r(:,j)=mat'*U(:,j)-beta(j)*V(:,j-1);
        for i = 1:j-1
            vr = i_hp_matrix_mul(complex(V(:,i)'),complex(r(:,j)));
            vrv = i_hp_matrix_mul_e(complex(V(:,i)),complex(vr));
            r(:,j) = i_hp_matrix_sub(complex(r(:,j)),complex(vrv));
%             r(:,j)=r(:,j)-(V(:,i)'*r(:,j))*V(:,i);
        end
        alpha(j) = nx(r(:,j));
        V(:,j) = i_hp_matrix_div_e(complex(r(:,j)),complex(alpha(j)));
%         V(:,j)=r(:,j)/alpha(j);
        mv = i_hp_matrix_mul(complex(mat),complex(V(:,j)));
        aU = i_hp_matrix_mul_e(complex(U(:,j)),complex(alpha(j)));
        P(:,j) = i_hp_matrix_sub(complex(mv),complex(aU));
%         P(:,j)=mat*V(:,j)-alpha(j)*U(:,j);
        for i = 1:j
            UP = i_hp_matrix_mul(complex(U(:,i)'),complex(P(:,j)));
            UPU = i_hp_matrix_mul_e(complex(U(:,i)),complex(UP));
            P(:,j) = i_hp_matrix_sub(complex(P(:,j)),complex(UPU));
%             P(:,j)=P(:,j)-(U(:,i)'*P(:,j))*U(:,i);
        end
        beta(j+1) = nx(P(:,j));
        if j~= 8
            U(:,j+1) = i_hp_matrix_div_e(complex(P(:,j)),complex(beta(j+1)));
        end
%         U(:,j+1)=P(:,j)/beta(j+1);
    else
        r(:,1) = i_hp_matrix_mul(complex(mat'),complex(U(:,1)));
        % mat_dec = mat;
%         r(:,1)=mat'*U(:,1);
        alpha(j) = nx(r(:,1));
        V(:,j) = i_hp_matrix_div_e(complex(r(:,1)),complex(alpha(j)));
%         V(:,j)=r(:,1)/alpha(j);
        mv = i_hp_matrix_mul(complex(mat),complex(V(:,j)));
        aU = i_hp_matrix_mul_e(complex(U(:,j)),complex(alpha(j)));
        P(:,j) = i_hp_matrix_sub(complex(mv),complex(aU));
%         P(:,j)=mat*V(:,j)-alpha(j)*U(:,j);
        UP = i_hp_matrix_mul(complex(U(:,j)'),complex(P(:,j)));
        UPU = i_hp_matrix_mul_e(complex(U(:,j)),complex(UP));
        P(:,j) = i_hp_matrix_sub(complex(P(:,j)),complex(UPU));
%         P(:,j)=P(:,j)-(U(:,j)'*P(:,j))*U(:,j);
        beta(j+1) = nx(P(:,j));
        U(:,j+1) = i_hp_matrix_div_e(complex(P(:,j)),complex(beta(j+1)));
%         U(:,j+1) = P(:,j)/beta(j+1);
    end
end
for i = 1:mini
    B(i,i) = alpha(i);
    B(i+1,i) = beta(i+1);
end
U = U(:,1:mini);
B = B(1:mini,:);
if isT==1
    T = U;
    B = B';
    U = V;
    V = T;
else
end
B = real(B);
end