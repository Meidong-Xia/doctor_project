function [U,BI,V]=bi_diag_svd(U,BI,V)
% BI = A;

% underflow = 1e-3;
tol = 1e-2;
epcl = 1e-6;
iter_num = 0;
% old_i_low = -1;
% old_i_high = -1;
[SIZEm,n] = size(BI);
fudge = min(SIZEm,n);
% maxit = 3*n*n;
maxit=30;
[low_bound_sigma,lad,~] = estimation(BI);
max_bound_sigma = max(abs(BI),[],'all');
% thresh = max(tol*low_bound_sigma,maxit*underflow);
while true
    iter_num = iter_num+1;
    if iter_num > maxit
        break
    end
    i_u = n-1;
    while((i_u >= 1)&&(abs(BI(i_u,i_u+1)) <= 1e-4))
        i_u = i_u-1;
    end
    if i_u == 0
        break
    end
    i_l = i_u-1;
    if i_l ~= 0
        while((abs(BI(i_l,i_l+1)) > 1e-4)&&(i_u >= 1))
            i_l = i_l-1;
            if i_l == 0
                break
            end
        end
    end
    if i_u == i_l+1
        % keyboard;
        [BI(i_l+1:i_u+1,i_l+1:i_u+1),U(:,i_l+1),U(:,i_u+1),V(:,i_l+1),V(:,i_u+1)] = mul_22_submatrix(BI(i_l+1:i_u+1,i_l+1:i_u+1),U(:,i_l+1),U(:,i_u+1),V(:,i_l+1),V(:,i_u+1));
        
        continue
    end
    % if(old_i_low~=i_l || old_i_high~=i_u)
    %     if(abs(BI(i_l+1,i_l+1)) >= abs(BI(i_u+1,i_u+1)))
    %         direction = 1;
    %     else
    %         direction = 2;
    %     end
    % end
    % old_i_low = i_l;
    % old_i_high = i_u;
    direction = 1;

    if direction == 1
        bdl = i_hp_div(abs(BI(i_u,i_u+1)),lad(i_u+1));
        if(bdl <= tol)
            BI(i_u,i_u+1) = 0;
        end
    end
    %compute shift
    ft = i_hp_mul(fudge,tol);
    ftl = i_hp_mul(ft,low_bound_sigma);
    ftldm = i_hp_div(ftl,max_bound_sigma);
%     fudge*tol*low_bound_sigma/max_bound_sigma
    if (ftldm <= epcl)
        shift = 0;
    else
        if direction == 1
            s = real(BI(i_u+1,i_u+1));
            r1 = real(i_hp_mul(BI(i_u,i_u),BI(i_u,i_u)));
            r2 = real(i_hp_mul(BI(i_u-1,i_u),BI(i_u-1,i_u)));
            r12 = real(i_hp_add(r1,r2));
            r3 = real(i_hp_mul(BI(i_u+1,i_u+1),BI(i_u+1,i_u+1)));
            r4 = real(i_hp_mul(+BI(i_u,i_u+1),BI(i_u,i_u+1)));
            r34 = real(i_hp_add(r3,r4));
            r12sr34 = real(i_hp_sub(r12,r34));
            d = real(i_hp_div(r12sr34,2));
%             d=((BI(i_u,i_u)*BI(i_u,i_u)+BI(i_u-1,i_u)*BI(i_u-1,i_u))-(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1)))/2;
            r1 = real(i_hp_mul(BI(i_u+1,i_u+1),BI(i_u+1,i_u+1)));
            r2 = real(i_hp_mul(BI(i_u,i_u+1),BI(i_u,i_u+1)));
            r12 = real(i_hp_add(r1,r2));
            r12d = real(i_hp_add(r12,d));
            d2 = real(i_hp_mul(d,d));
            r3 = real(i_hp_mul(BI(i_u,i_u),BI(i_u,i_u)));
            r4 = real(i_hp_mul(r3,BI(i_u,i_u+1)));
            r5 = real(i_hp_mul(r4,BI(i_u,i_u+1)));
            r6 = real(i_hp_add(d2,r5));
            sqr6 = real(i_hp_sqrt(r6));
            ssqr6 = real(i_hp_mul(sign(d),sqr6));
            shift = real(i_hp_sub(r12d,ssqr6));
%             shift=(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1))+d-sign(d)*sqrt(d*d+BI(i_u,i_u)*BI(i_u,i_u)*BI(i_u,i_u+1)*BI(i_u,i_u+1));
        else
            s = real(BI(i_l+1,i_l+1));
            r1 = real(i_hp_mul(BI(i_l+2,i_l+2),BI(i_l+2,i_l+2)));
            r2 = real(i_hp_mul(BI(i_l+2,i_l+3),BI(i_l+2,i_l+3)));
            r12 = real(i_hp_add(r1,r2));
            r3 = real(i_hp_mul(BI(i_l+1,i_l+1),BI(i_l+1,i_l+1)));
            r4 = real(i_hp_mul(BI(i_l+1,i_l+2),BI(i_l+1,i_l+2)));
            r34 = real(i_hp_add(r3,r4));
            r12sr34 = real(i_hp_sub(r12,r34));
            d = real(i_hp_div(r12sr34,2));
%             d=((BI(i_l+2,i_l+2)^2+BI(i_l+2,i_l+3)^2)-(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2))/2;
            r1 = real(i_hp_mul(BI(i_l+1,i_l+1),BI(i_l+1,i_l+1)));
            r2 = real(i_hp_mul(BI(i_l+1,i_l+2),BI(i_l+1,i_l+2)));
            r12 = real(i_hp_add(r1,r2));
            r12d = real(i_hp_add(r12,d));
            d2 = real(i_hp_mul(d,d));
            r3 = real(i_hp_mul(BI(i_l+2,i_l+2),BI(i_l+2,i_l+2)));
            r4 = real(i_hp_mul(BI(i_l+1,i_l+2),BI(i_l+1,i_l+2)));
            r5 = real(i_hp_mul(r3,r4));
            dr5 = real(i_hp_add(d2,r5));
            sqdr5 = real(i_hp_sqrt(dr5));
            ssqdr5 = real(i_hp_mul(sign(d),sqdr5));
            shift = real(i_hp_sub(r12d,ssqdr5));
%             shift=(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2)+d-sign(d)*sqrt(d^2+BI(i_l+2,i_l+2)^2*BI(i_l+1,i_l+2)^2);
        end
        shift2 = real(i_hp_mul(shift,shift));
        s2 = real(i_hp_mul(s,s));
        shiftds = real(i_hp_div(shift2,s2));
        if(shiftds)<=epcl
            shift = 0;
        end
    end
    if shift ~= 0
        if direction == 1
            % BIO=BI;
            % UO=U;
            % VO=V;
            % keyboard;
            [BI,U,V] = QR_Wilkinson_shift_Iteration_once(BI,U,V,i_l,i_u,shift);
            % norm(U*BI*V'-UO*BIO*VO',"fro")/norm(UO*BIO*VO',"fro");
        else
            [BI,U]=QR_Wilkinson_shift_Iteration_once_upward(BI,U,i_l,i_u,shift);
        end
    else
        if direction==1
            [BI,U,V]=QR_zero_shift_once_iteration(BI,U,V,i_l,i_u);
        else
            [BI,U]=QR_zero_shift_once_iteration_upward(BI,U,i_l,i_u);
        end
    end
end
V=V';
end
