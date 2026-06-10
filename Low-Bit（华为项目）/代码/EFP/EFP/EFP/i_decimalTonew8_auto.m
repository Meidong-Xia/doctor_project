function [fp8Binary,result_config] = i_decimalTonew8_auto(decimalValue,config,base,fraction_tables)
% config = [e_bit m_bit_real bias_real e_bit_image m_bit_image bias_image] 
    % base = 10;
    % config = {[2 5 2 2 5 2]} ;
    % decimalValue = 0.000124+234532i;

    % global mu_xishu;
    % global m;
    config_double = config{1};
    dec_real = real(decimalValue);
    dec_imag = imag(decimalValue);
    [fp8Binary_real,result_config_real] = decimalTonew8_auto(dec_real,config_double(1:3),base,fraction_tables);
    [fp8Binary_imag,result_config_imag] = decimalTonew8_auto(dec_imag,config_double(6:8),base,fraction_tables);

    % fp8_real = new8Todecimal(fp8Binary_real,result_config_real);
    % fp8_imag = new8Todecimal(fp8Binary_imag,result_config_imag);
    % config_double(4) = mu_xishu*((fp8_real - dec_real) / dec_real);
    % config_double(9) = mu_xishu*((fp8_imag - dec_imag) / dec_imag);
    % config_double(4) = mu_xishu*(rand-0.5)*2^(-m);
    % config_double(9) = mu_xishu*(rand-0.5)*2^(-m);
    
    fp8Binary = [fp8Binary_real fp8Binary_imag];

    result_config = [result_config_real config_double(4:5) result_config_imag config_double(9:10)];

end