function result = fp8Addition_e5m2(fp8num1, fp8num2)
    % 输入参数 fp8num1 和 fp8num2 为两个 8 位的 fp8 数字的二进制字符串
    
    % 解析输入，分别提取符号位、指数位和小数位
    sign1 = str2double(fp8num1(1));
    sign2 = str2double(fp8num2(1));
    
    exponent1 = bin2dec(fp8num1(2:6));
    exponent2 = bin2dec(fp8num2(2:6));
    
    fraction1 = bin2dec(fp8num1(7:8));
    fraction2 = bin2dec(fp8num2(7:8));
    
    % 对齐指数位
    if exponent1 > exponent2
        fraction2 = fraction2 * 2^(exponent1 - exponent2);
        exponent2 = exponent1;
    else
        fraction1 = fraction1 * 2^(exponent2 - exponent1);
        exponent1 = exponent2;
    end
    
    % 根据符号位调整小数部分
    if sign1 == 1
        fraction1 = -fraction1;
    end
    if sign2 == 1
        fraction2 = -fraction2;
    end
    
    % 小数位相加
    resultFraction = fraction1 + fraction2;
    
    % 判断进位
    if resultFraction > 3
        resultFraction = resultFraction - 4;
        exponent1 = exponent1 + 1;
    elseif resultFraction < -4
        % 负溢出
        resultSign = 1;
        resultExponent = 31; % 最大指数
        resultFraction = 0;  % 最小小数
    end
    if resultFraction >= 4
        resultFraction = resultFraction - 4;
        exponent1 = exponent1 + 1;
    end
    
    % 判断溢出
    if exponent1 > 30
        result = 'Overflow';
        return;
    end
    
    % 计算新的指数位
    resultExponent = exponent1;
    
    % 判断符号位
    if resultFraction < 0
        resultSign = 1;
        resultFraction = -resultFraction;
    else
        resultSign = 0;
    end
    
    % 转换小数位为二进制字符串
    resultFractionBinary = dec2bin(resultFraction, 2);
    
    % 补齐小数位长度
    while length(resultFractionBinary) < 2
        resultFractionBinary = ['0', resultFractionBinary];
    end
    
    % 拼接结果
    result = [num2str(resultSign), dec2bin(resultExponent, 5), resultFractionBinary];
end
