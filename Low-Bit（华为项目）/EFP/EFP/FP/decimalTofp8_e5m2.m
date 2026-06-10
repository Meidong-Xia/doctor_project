function fp8Binary = decimalTofp8_e5m2(decimalValue)
%         decimalValue = 2^-16;
    bias = 15;
    % 确保指数在合法范围内,如果大于等于0 11110 111，则溢出
    if decimalValue >= 61440 
        fp8Binary = '01111100';
        return;
    elseif decimalValue <= -61440
        fp8Binary = '11111100';
        return;
    elseif isnan(decimalValue)
        % NaN
        fp8Binary = '01111101'; 
        return;
    % 处理特殊值
    elseif isinf(decimalValue)
        if decimalValue < 0
            % 负无穷
            fp8Binary = '11111100';
        else
            % 正无穷
            fp8Binary = '01111100';
        end
        return;
    end
        
    hex = '0123456789abcdef'; %Hex characters
    h = num2hex(decimalValue);	%Convert from float to hex characters
    hc = num2cell(h); %Convert to cell array of chars
    nums =  cellfun(@(x) find(hex == x) - 1, hc); %Convert to array of numbers
    bins = dec2bin(nums, 4); %Convert to array of binary number strings
    fp64Binary = reshape(bins.', 1, numel(bins)); %Reshape into horizontal vector
    % 输入参数 fp64Binary 应为一个包含fp64的二进制字符串
    
    % 提取符号、指数和小数位
    signBit = fp64Binary(1);
    exponentBits = fp64Binary(2:12);
    fractionBits = fp64Binary(13:end);

    % 将fp64的指数部分减去1023，以得到fp8的指数部分
    exponent = bin2dec(exponentBits) - 1023;

    % 处理规范数和非规范数
    if exponent < -14  % 非规范数
        % 将fp64小数部分扩大到 2^2 的范围，使其变为fp8规范数
        fraction = abs(decimalValue) / 2^(-14);
        % 添加额外的条件判断
        if fraction >= 0.875
            % 大于等于111的情况
            fp8Binary = [signBit '0000100'] ;
        elseif fraction > 0.625
            % 大于等于0.625的情况
            fp8Binary = [signBit '0000011'] ;
        elseif fraction >= 0.375
            % 大于等于0.375的情况
            fp8Binary = [signBit '0000010'] ;
        elseif fraction > 0.125
            % 大于等于0.125的情况
            fp8Binary = [signBit '0000001'] ;
        else
            % 小于0.125的情况
            fp8Binary = [signBit '0000000'] ;
        end
    else  % 规范数
        % 提取小数部分的前两位（四舍五入）
        fraction = fractionBits(1:2);
        twoBit = fractionBits(2);
        threeBit = fractionBits(3);
        % 判断 fractionBits 的第4位到最后一位是否存在1
        hasOnes = any(fractionBits(4:end) == '1');
        quan1 = '11';
        if threeBit == '1'
            if isequal(fractionBits(1:2), quan1)
                fraction = '00';
                exponent = exponent + 1;
            elseif hasOnes == 1 || twoBit == '1'
                fraction = dec2bin(bin2dec(fraction) + 1, 2);
            end
        end
        % 构建fp8二进制表示
        fp8Binary = [signBit, dec2bin(exponent+bias, 5), fraction];
    end
end


