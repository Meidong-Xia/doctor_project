function [eta,F] = f_mmse_precoding(H,P,Nu,Ns,Nr,snr,flg_b);
if flg_b == 1
    F = H'*inv(H*H'+1/snr/P*eye(Nr));
    mask = eye(Nr);
    mask((Nu*Ns+1):Nr,(Nu*Ns+1):Nr) = 0;
    eta = sqrt(P/trace(F*mask*F'));
    F = eta*F;
elseif flag_b == 2
    [~,~,tmp] = i_hp_div(1,snr);
    [~,~,coeff_value] = i_hp_div(tmp,P);
    coeff_matrix = coeff_value*eye(Nr);
    [~,~,tmp1] = i_hp_matrix_mul(complex(H),complex(H'));
    [~,~,tmp2] = i_hp_matrix_add(complex(coeff_matrix),complex(tmp1));
    [~,~,F] = i_hp_matrix_mul(complex(H'),complex(i_GaussianElimination_half( tmp2)));
    mask = eye(Nr);
    mask((Nu*Ns+1):Nr,(Nu*Ns+1):Nr) = 0;
    [~,~,F_tmp] = i_hp_matrix_mul(complex(mask),complex(F'));
    [~,~,F_head] = i_hp_matrix_mul(complex(F),complex(F_tmp));
    power = 0;
    for k=1:size(F_head,1)
        [~,~,power] = i_hp_add(power,F_head(k,k));
    end
    [~,~,tmp3] = i_hp_div(P,power);
    [~,eta] = r_hp_sqrt(tmp3);
    [~,~,F] = i_hp_matrix_mul(complex(eta*diag(ones(size(F,1),1))),complex(F));
end

end