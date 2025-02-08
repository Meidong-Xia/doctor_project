function decimalValue = new8Todecimal(fp8Binary,config, base, fraction_tables)
    % 输入参数 fp8Binary 应为一个包含8位二进制字符串的字符数组
     % bias = 2^(e_bit-1);
     %  fp8Binary = '010000001';
     %  e_bit = 3;
     %  m_bit = 5;
     %  base = 2;
    % 验证输入长度
    % global base;
    e_bit = config(1);
    m_bit = config(2);
    bias = config(3);
    if numel(fp8Binary) ~= e_bit+m_bit+1
        error('输入的二进制字符串不符。');
    end
     % fp8Binary = '01011111';
    % 提取符号、指数和小数位
    signBit = str2double(fp8Binary(1));
    exponentBits = bin2dec(fp8Binary(2:e_bit+1));
    fractionBits = bin2dec(fp8Binary(e_bit+2:e_bit+1+m_bit));
    
    if strcmp(fp8Binary,['1',dec2bin(0,m_bit+e_bit)])
        decimalValue = NaN;
        return;
    elseif strcmp(fp8Binary,dec2bin(0,m_bit+e_bit+1)) 
        decimalValue = 0;
        return;
    end

    % bias = 2;
    exponent = exponentBits - bias;  

    % 计算尾数的实际值
    % fraction_table = [0 1 15/14 15/13 5/4 100/73 3/2 8/5 17/10 9/5 19/10 2 20/9 40/17 2.5 14/5 10^(1/2) 50/14 4 17/4 9/2 5 100/19 50/9 100/17 6.25 20/3 73/10 8 26/3 28/3 Inf];
    % fraction_table = [0 1 10/9 5/4 100/73 3/2 8/5 17/10 9/5 19/10 2 40/19 20/9 40/17 2.5 14/5 10^(1/2) 50/14 4 17/4 9/2 19/4 5 100/19 50/9 100/17 6.25 20/3 73/10 8 9 Inf];
    
    % fraction_table = zeros(1,2^m_bit);
    % base = 2;
    % parfor k = 1:2^m_bit
    %     fraction_table(k) = base^((k-1)/2^m_bit);
    % end

    % global fraction_tables;
    fraction_tables = fraction_tables.Value;
    if base == 10
        fraction_table = fraction_tables{2, m_bit};
    else
        fraction_table = fraction_tables{1, m_bit};
    end



    fraction = fraction_table(fractionBits+1);

    % 计算浮点数的十进制值

    decimalValue = (-1)^signBit * base^exponent * fraction;
    
end
