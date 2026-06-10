function [fp8Binary,result_config] = decimalTonew8_auto(decimalValue,config,base,fraction_tables)
    % 提取底数和指数位
    % decimalValue = -0.0219;
    % config = [2 5 2];
    % global base;
    % e_bit = 3;
    %  bias = 2^(e_bit-1);
     
     
     % m_bit = 6;
    % fraction_table = [0 1 15/14 15/13 5/4 100/73 3/2 8/5 17/10 9/5 19/10 2 20/9 40/17 2.5 14/5 10^(1/2) 50/14 4 17/4 9/2 5 100/19 50/9 100/17 6.25 20/3 73/10 8 26/3 28/3 Inf];
    % fraction_table = [0 1 10/9 5/4 100/73 3/2 8/5 17/10 9/5 19/10 2 40/19 20/9 40/17 2.5 14/5 10^(1/2) 50/14 4 17/4 9/2 19/4 5 100/19 50/9 100/17 6.25 20/3 73/10 8 9 Inf];
    


    e_bit = config(1);
    m_bit = config(2);
    bias = config(3);

    result_config(1) = config(1);
    result_config(2) = config(2);
    result_config(3) = config(3);

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

    if decimalValue < 0
        signBit = '1';
    else 
        signBit = '0';
    end

   if isnan(decimalValue)||isinf(decimalValue)
       fp8Binary = ['1',dec2bin(0,e_bit+m_bit)];%
       return ;
   elseif decimalValue == 0
       fp8Binary = dec2bin(0,e_bit+m_bit+1);
       return ;
   end

    mantissa = abs(decimalValue) / base^floor(log(abs(decimalValue)) / log(base));

    exponent = floor(log(abs(decimalValue)) / log(base));

    for i = 1:2^m_bit
        if i == 2^m_bit
            if mantissa < (fraction_table(i)+base)/2
                mantissa_bit = dec2bin(i-1,m_bit);
            else 
                mantissa_bit = dec2bin(0,m_bit);
                exponent = exponent+1;
            end
        elseif mantissa < (fraction_table(i)+fraction_table(i+1))/2
            mantissa_bit = dec2bin(i-1,m_bit);
            break
        end
    end

    exponent = exponent + bias;
   
    if exponent > 2^e_bit-1
        result_config(3) = result_config(3) - (exponent-(2^e_bit-1));
        exponent = 2^e_bit-1;
    elseif exponent < 0
        result_config(3) = result_config(3) - (exponent);
        exponent = 0;
    end

    if exponent == 0 && bin2dec(mantissa_bit) == 0
           exponent = 1;
           result_config(3) = result_config(3) +1;
    end

    fp8Binary = [signBit, dec2bin(exponent, e_bit), mantissa_bit];
end


