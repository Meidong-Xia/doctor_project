function [result_fp8Binary, result_config] = EFP_mul(fp8Binary1, config1, fp8Binary2,config2)
    % % 检查输入的长度是否正确
    % fp8Binary1 = '00000001';
    % fp8Binary2 ='00011111';
    % e_bit1 = 2;
    % e_bit2 = 2;
    % m_bit1 = 5; 
    % m_bit2 = 5; 
    % bias1 = 2;
    % bias2 = 2;
    % base = 2;

    e_bit1 = config1(1);
    m_bit1 = config1(2);
    bias1 = config1(3);
    e_bit2 = config2(1);
    m_bit2 = config2(2);
    bias2 = config2(3);

    result_config = zeros(1,3);
    result_config(1) = e_bit1;
    result_config(2) = m_bit1;
    result_m_bit = m_bit1;
    if numel(fp8Binary1) ~= e_bit1 + m_bit1 + 1 || numel(fp8Binary2) ~= e_bit2 + m_bit2 + 1
        keyboard;
        error('输入的二进制字符串长度不正确。');
    end
    
    % 判断是否为NaN或者乘0
    if strcmp(fp8Binary1,['1',dec2bin(0,m_bit1+e_bit1)]) || strcmp(fp8Binary2,['1',dec2bin(0,m_bit2+e_bit2)])
        result_fp8Binary = ['1',dec2bin(0,result_m_bit+e_bit1)];
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    elseif strcmp(fp8Binary1,dec2bin(0,m_bit1+e_bit1+1)) || strcmp(fp8Binary2,dec2bin(0,m_bit2+e_bit2+1))
        result_fp8Binary = dec2bin(0,1+e_bit1+result_m_bit);
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    end

    % 提取符号、指数和尾数
    signBit1 = str2double(fp8Binary1(1));
    exponentBits1 = bin2dec(fp8Binary1(2:e_bit1+1));
    fractionBits1 = bin2dec(fp8Binary1(e_bit1+2:e_bit1+m_bit1+1));
    
    signBit2 = str2double(fp8Binary2(1));
    exponentBits2 = bin2dec(fp8Binary2(2:e_bit2+1));
    fractionBits2 = bin2dec(fp8Binary2(e_bit2+2:e_bit2+m_bit2+1));

    % 确保尾数位宽度相同
    if m_bit1 < m_bit2
        % 如果fp8Binary1的尾数位较短，补0至与fp8Binary2相同的位宽
        fractionBits1 = fractionBits1 * 2^(m_bit2 - m_bit1);
    elseif m_bit1 > m_bit2
        % 如果fp8Binary2的尾数位较短，补0至与fp8Binary1相同的位宽
        fractionBits2 = fractionBits2 * 2^(m_bit1 - m_bit2);
    end
    
    
    % 计算指数和尾数的值
    exponent1 = exponentBits1 - bias1;
    exponent2 = exponentBits2 - bias2;
    

    % 计算乘积的指数和尾数
    exponent_sum = exponent1 + exponent2;
    fraction_product = fractionBits1 + fractionBits2;
    
    % 将乘积尾数转换为二进制，尾数位宽选择最大尾数位宽
    max_m_bit = max(m_bit1, m_bit2);
    fraction_product_binary = dec2bin(fraction_product, max_m_bit);
    
    % 如果尾数位宽超过两个加数的最大尾数位宽，则指数位加1，尾数舍去最高位
    if length(fraction_product_binary) > max_m_bit
        fraction_product_binary = fraction_product_binary(2:end);
        exponent_sum = exponent_sum + 1;
    end

    result_config(3) = 2;
    exponent_sum = exponent_sum + result_config(3);

    % 判断是否需要调整指数和尾数
    if exponent_sum > (2^result_config(1) - 1)
        result_config(3) = result_config(3) - (exponent_sum - (2^result_config(1) - 1));
        exponent_sum = (2^result_config(1) - 1);
        
    elseif exponent_sum < 0
        result_config(3) = result_config(3) - (exponent_sum);
        exponent_sum = 0;
    end
    
    if exponent_sum == 0 && bin2dec(fraction_product_binary) == 0
        exponent_sum = exponent_sum+1;
        result_config(3) = result_config(3)+1;
    end

    % 更新尾数位宽
    fraction_binary = fraction_product_binary;

    
    
    % 计算符号位
    sign_result = xor(signBit1, signBit2);

    if sign_result == 1 &&  bin2dec(fraction_binary) == 0  && exponent_sum == 0
        sign_result = 0 ;
    end

    % 返回乘积结果
    result_fp8Binary = [num2str(sign_result), dec2bin(exponent_sum, result_config(1)), fraction_binary];

end
