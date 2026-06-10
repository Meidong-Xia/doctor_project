function Plist = f_power_allocation_1(U,Sigma,V,Ns,P,m_bit, base, fraction_tables,table)

V_tmp = EFP_mul_matrix_2(V, V', m_bit, base, fraction_tables,table);
P_tmp = 0;
for i = 1:size(V_tmp,1)
    P_tmp = EFP_add_matrix_2(P_tmp,V_tmp(i,i), m_bit, base, fraction_tables,table);
end
P_sqrt = EFP_div_matrix_e_2(P,P_tmp,m_bit, base, fraction_tables,table);
P_tmp = i_EFP_sqrt_2(P_sqrt, m_bit, base, fraction_tables,table);
Plist = P_tmp*single(diag(ones(Ns,1)));

% [V_efp, V_config] = decToEFP_auto(V,m_bit,base, fraction_tables);
% [V_trans_efp, V_trans_config] = hermitian_transpose(V_efp,V_config);
% 
% [V_tmp_efp, V_tmp_config] = EFP_mul_matrix(V_efp,V_config,V_trans_efp,V_trans_config, base, ...
%     fraction_tables,table);
% 
% [P_tmp_efp,P_tmp_efp_config] = decToEFP_auto(0,m_bit,base, fraction_tables);
% [P_efp,P_config] = decToEFP_auto(P,m_bit,base, fraction_tables);
% 
% for i = 1:size(V,1)
%     [P_tmp_efp, P_tmp_efp_config] = EFP_add_matrix(P_tmp_efp,P_tmp_efp_config,V_tmp_efp(i,i), ...
%         V_tmp_config(i,i), ...
%         base, fraction_tables,table);
% end
% 
% [P_sqrt_efp, P_sqrt_efp_config] = EFP_div_matrix_e(P_efp,P_config,P_tmp_efp,P_tmp_efp_config,base, ...
%     fraction_tables,table);
% 
% [P_tmp,P_tmp_config] = i_EFP_sqrt(P_sqrt_efp,P_sqrt_efp_config, base, ...
%     fraction_tables,table);
% 
% Plist = EFPTodec({P_tmp}, {P_tmp_config},base,fraction_tables);
% Plist = Plist*(diag(ones(Ns,1)));

end