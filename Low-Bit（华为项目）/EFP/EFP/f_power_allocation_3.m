function Plist = f_power_allocation_3(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho)
    Sigma = EFP_mul_matrix_2(Sigma,Sigma,m_bit, base, fraction_tables,table);
    b = 0;
    for i = 1:Ns
        tmp = EFP_div_matrix_e_2(sigma_n,Sigma(i,i),m_bit, base, fraction_tables,table);
        b = EFP_add_matrix_2(b,tmp,m_bit, base, fraction_tables,table);
    end
    tmp2 = EFP_mul_matrix_2(Ns,P,m_bit, base, fraction_tables,table);
    tmp2 = EFP_mul_matrix_2(tmp2,rho,m_bit, base, fraction_tables,table);

    tmp3 = EFP_add_matrix_2(b,tmp2,m_bit, base, fraction_tables,table);

    tmp4 = EFP_mul_matrix_2(P,P,m_bit, base, fraction_tables,table);

    mu = EFP_div_matrix_e_2(tmp3,tmp4,m_bit, base, fraction_tables,table);

    Plist = zeros(Ns,1);

    for i = 1:Ns
        tmp1 = EFP_mul_matrix_2(P,rho,m_bit, base, fraction_tables,table);
        tmp2 = EFP_div_matrix_e_2(sigma_n,Sigma(i,i),m_bit, base, fraction_tables,table);
        tmp3 = EFP_add_matrix_2(tmp1,tmp2,m_bit, base, fraction_tables,table);
        tmp4 = EFP_mul_matrix_2(mu,P,m_bit, base, fraction_tables,table);
        Plist(i,1) = EFP_div_matrix_e_2(tmp3,tmp4,m_bit, base, fraction_tables,table);
    end
    Plist = diag(Plist);
end