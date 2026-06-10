function [eta,F] = f_zf_precoding(H,P,Nu,Ns,Nr,flg_b)
if flg_b == 1 % 64bit
    F = H'*inv(H*H');
    mask = eye(Nu*Nr);
    for i=1:Nu
        mask(((i-1)*Nr+Ns+1):(i*Nr),((i-1)*Nr+Ns+1):(i*Nr)) = 0;
    end
    eta = sqrt(P/trace(F*mask*F'));
    F = eta*F;
elseif flg_b == 2 % 32bit
    F = single(H')*GaussianElimination_single(single(H)*single(H'));
    mask = eye(Nu*Nr);
    for i=1:Nu
        mask(((i-1)*Nr+Ns+1):(i*Nr),((i-1)*Nr+Ns+1):(i*Nr)) = 0;
    end
    eta = sqrt(single(P)/single(trace(single(F)*single(mask)*single(F'))));
    F = single(eta)*single(F);
elseif flg_b == 3 % 16bit
    [~,~,tmp1] = i_hp_matrix_mul(complex(H),complex(H'));
    [~,~,F] = i_hp_matrix_mul(complex(H'),complex(improve_GaussianElimination_half(tmp1)));
    mask = eye((Nu*Nr));
    for i=1:Nu
        mask(((i-1)*Nr+Ns+1):(i*Nr),((i-1)*Nr+Ns+1):(i*Nr)) = 0;
    end
    [~,~,F_tmp] = i_hp_matrix_mul(complex(mask),complex(F'));
    [~,~,F_head] = i_hp_matrix_mul(complex(F),complex(F_tmp));
    power = 0;
    for k=1:size(F_head,1)
        [~,~,power] = i_hp_add(power,F_head(k,k));
    end
    [~,~,tmp2] = i_hp_div(P,power);
    [~,eta] = i_hp_sqrt(tmp2);
    [~,~,F] = i_hp_matrix_mul(complex(eta*diag(ones(size(F,1),1))),complex(F));
end
end