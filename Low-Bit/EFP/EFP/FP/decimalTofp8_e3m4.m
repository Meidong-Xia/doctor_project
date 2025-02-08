function fp8Binary = decimalTofp8_e3m4(decimalValue,bias)
%      decimalValue = (1+1/8+1/16)*2^-5;
    % bias = 3;
    % 确保指数在合法范围内,如果大于等于0 11110 111，则溢出
    if decimalValue >= 2^(7-bias)*(1+1/2+1/4+1/8)
        fp8Binary = '01111110';
        return;
    elseif decimalValue <= -2^(7-bias)*(1+1/2+1/4+1/8)
        fp8Binary = '11111110';
        return;
    elseif isnan(decimalValue)
        % NaN
        fp8Binary = '01111111'; 
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
    if exponent < 1-bias  % 非规范数
        % 将fp64小数部分扩大到 2^2 的范围，使其变为fp8规范数
        fraction = abs(decimalValue) / 2^(1-bias);
        % 添加额外的条件判断
        if fraction >= (31/32)
            % 大于等于1111的情况
            fp8Binary = [signBit '0010000'] ;
        elseif fraction > 29/32
            % 大于等于1101的情况
            fp8Binary = [signBit '0001111'] ;
        elseif fraction >= 27/32
            % 大于等于1011的情况
            fp8Binary = [signBit '0001110'] ;
        elseif fraction > 25/32
            % 大于等于1001的情况
            fp8Binary = [signBit '0001101'] ;
        elseif fraction >= 23/32
            % 大于等于0111的情况
            fp8Binary = [signBit '0001100'] ;
        elseif fraction > 21/32
            % 大于等于0101的情况
            fp8Binary = [signBit '0001011'] ;
        elseif fraction >= 19/32
            % 大于等于0011的情况
            fp8Binary = [signBit '0001010'] ;
        elseif fraction > 17/32
            % 大于等于0001的情况
            fp8Binary = [signBit '0001001'] ;
        elseif fraction > 15/32
            % 大于等于1101的情况
            fp8Binary = [signBit '0001000'] ;
        elseif fraction > 13/32
            % 大于等于1101的情况
            fp8Binary = [signBit '0000111'] ;
        elseif fraction >= 11/32
            % 大于等于1011的情况
            fp8Binary = [signBit '0000110'] ;
        elseif fraction > 9/32
            % 大于等于1001的情况
            fp8Binary = [signBit '0000101'] ;
        elseif fraction >= 7/32
            % 大于等于0111的情况
            fp8Binary = [signBit '0000100'] ;
        elseif fraction > 5/32
            % 大于等于0101的情况
            fp8Binary = [signBit '0000011'] ;
        elseif fraction >= 3/32
            % 大于等于0011的情况
            fp8Binary = [signBit '0000010'] ;
        elseif fraction > 1/32
            % 大于等于0001的情况
            fp8Binary = [signBit '0000001'] ;
        else
            % 小于0.125的情况
            fp8Binary = [signBit '0000000'] ;
        end
    else  % 规范数
        % 提取小数部分的前两位（四舍五入）
        fraction = fractionBits(1:4);
        threeBit = fractionBits(4);
        fourBit = fractionBits(5);
        % 判断 fractionBits 的第5位到最后一位是否存在1
        hasOnes = any(fractionBits(6:end) == '1');
        quan1 = '1111';
        if fourBit == '1'
            if isequal(fractionBits(1:4), quan1)
                fraction = '0000';
                exponent = exponent + 1;
            elseif hasOnes == 1 || threeBit == '1'
                fraction = dec2bin(bin2dec(fraction) + 1, 4);
            end
        end
        % 构建fp8二进制表示
        fp8Binary = [signBit, dec2bin(exponent+bias, 3), fraction];
    end
end


