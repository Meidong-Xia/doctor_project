function decimalValue = fp8Todecimal_e5m2(fp8Binary)
    % 输入参数 fp8Binary 应为一个包含8位二进制字符串的字符数组
    
    % 验证输入长度
    if numel(fp8Binary) ~= 8
        error('输入的二进制字符串必须为8位。');
    end
    % 提取符号、指数和小数位
    signBit = str2double(fp8Binary(1));
    exponentBits = bin2dec(fp8Binary(2:6));
    fractionBits = bin2dec(fp8Binary(7:8));

    % 判断特殊值
    if exponentBits == 31
        if fractionBits == 0
            % INF
            decimalValue = inf * (-1)^signBit;
            return;
        elseif fractionBits == 1 || fractionBits == 2 || fractionBits == 3
            % NaN
            decimalValue = NaN;
            return;
        end
    end

    % 计算指数和尾数的实际值
    exponent = exponentBits - 15;  % 偏移为15
    fraction = fractionBits / 2^2;  % 因为有2位小数

    % 计算浮点数的十进制值
    if exponent == -15
        % subnormal
        decimalValue = (-1)^signBit * 2^(-14) * fraction;
    else
        % normal
        decimalValue = (-1)^signBit * 2^exponent * (1 + fraction);
    end
end
