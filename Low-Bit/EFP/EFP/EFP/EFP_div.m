function [result_fp8Binary, result_config] = EFP_div(fp8Binary1, config1, fp8Binary2, config2)
    % % 检查输入的长度是否正确
    %         fp8Binary1 = '00000001';
    %         fp8Binary2 = '11100000';
    % e_bit1 = 2;
    % e_bit2 = 2;
    % m_bit1 = 5; 
    % m_bit2 = 5; 
    % bias1 = 2;
    % bias2 = 2;
    % base = 10;

    e_bit1 = config1(1);
    bias1 = config1(3);
    e_bit2 = config2(1);
    m_bit2 = config2(2);

    result_config = zeros(1,3);
    result_config(1) = e_bit1;
    result_config(2) = config1(2);
    result_m_bit = m_bit2;
    
    
    % 提取符号、指数和尾数
    exponentBits2 = bin2dec(fp8Binary2(2:e_bit2+1));
    fractionBits2 = bin2dec(fp8Binary2(e_bit2+2:e_bit2+m_bit2+1));


    % 判断是否为NaN或者乘0
    % 判断是否为NaN或者乘0
    if  strcmp(fp8Binary2,['1',dec2bin(0,m_bit2+e_bit2)])||strcmp(fp8Binary2,dec2bin(0,m_bit2+e_bit2+1))
        result_fp8Binary = ['1',dec2bin(0,result_m_bit+e_bit1)];
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    end
    
    exp2 = dec2bin(3-exponentBits2,e_bit2);
    config2(3) = -(config2(3)-4);

    if fractionBits2 == 0 
        config2(3) = config2(3)-1;
        if bin2dec(exp2) == 0 
            exp2 = dec2bin(3-exponentBits2+1,e_bit2);
            config2(3) = config2(3)+1;
        end
    end

    m_temp = dec2bin(2^m_bit2-fractionBits2,m_bit2);
    m2 = m_temp(end-m_bit2+1:end);

    fp8Binary2(2:e_bit2+1) = exp2;
    fp8Binary2(e_bit2+2:e_bit2+m_bit2+1) = m2;

    [result_fp8Binary, result_config] = EFP_mul(fp8Binary1, config1, fp8Binary2, config2);
    % 返回乘积结果
    return;
end
