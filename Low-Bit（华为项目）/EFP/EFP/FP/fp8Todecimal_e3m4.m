function decimalValue = fp8Todecimal_e3m4(fp8Binary,bias)
    % 输入参数 fp8Binary 应为一个包含8位二进制字符串的字符数组
    % bias = 3;
    % 验证输入长度
    if numel(fp8Binary) ~= 8
        error('输入的二进制字符串必须为8位。');
    end
    
    % 提取符号、指数和小数位
    signBit = str2double(fp8Binary(1));
    exponentBits = bin2dec(fp8Binary(2:4));
    fractionBits = bin2dec(fp8Binary(5:8));

    % 判断特殊值
    if exponentBits == 7
        if fractionBits == 15
            % NaN
            decimalValue = NaN;
            return;
        else
            exponent = exponentBits - bias;  % 偏移为3
        end
    elseif exponentBits == 0
        if fractionBits == 0
            % 0
            decimalValue = 0;
            return;
        else
            % 非常规数
            exponent = -bias;  % subnormal
        end
    else
        % normal
        exponent = exponentBits - bias;  % 偏移为3
    end

    % 计算尾数的实际值
    fraction = fractionBits / 2^4;  % 因为有4位小数

    % 计算浮点数的十进制值
    if exponent == -bias
        % subnormal
        decimalValue = (-1)^signBit * 2^(1-bias) * fraction;
    else
        % normal
        decimalValue = (-1)^signBit * 2^exponent * (1 + fraction);
    end
end
