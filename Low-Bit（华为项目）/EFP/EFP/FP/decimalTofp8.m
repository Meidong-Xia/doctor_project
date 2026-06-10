function fp8Binary = decimalTofp8(decimalValue,config)
     %  decimalValue = 0;
     e_bit = config(1);
     m = config(2);
     bias = 2^(e_bit-1)-1;
     
     % m = 3;
    % 确保指数在合法范围内,如果大于等于0 11110 111，则溢出
    if decimalValue >= 2^(2^e_bit-1-bias)*(1+(2^(m-1)-1)/2^(m-1))
        fp8Binary = ['0',repmat('1',1,e_bit-1+m),'0'];
        return;
    elseif decimalValue <= -2^(2^e_bit-1-bias)*(1+(2^(m-1)-1)/2^(m-1))
        fp8Binary = ['1',repmat('1',1,e_bit-1+m),'0'];
        return;
    elseif isnan(decimalValue)
        % NaN
        fp8Binary = ['0',repmat('1',1,e_bit+m)]; 
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
        for i = flip(0:2^m-1)
            if mod(i, 2) == 1
                if fraction >= (i/2^m+(i+1)/2^m)/2
                    if i == 2^m-1
                        fp8Binary = [signBit, dec2bin(0,e_bit-1),'1',dec2bin(0,m)] ;
                        break;
                    else
                        fp8Binary = [signBit, dec2bin(0,e_bit),dec2bin(i+1,m)] ;
                        break;
                    end
                end
            else 
                if i == 0
                    if fraction > (i/2^m+(i+1)/2^m)/2
                        fp8Binary = [signBit, dec2bin(0,e_bit),dec2bin(i+1,m)] ;
                        break;
                    else
                        fp8Binary = [signBit, dec2bin(0,e_bit),dec2bin(i,m)] ;
                        break;
                    end
                else 
                    if fraction > (i/2^m+(i+1)/2^m)/2
                        fp8Binary = [signBit, dec2bin(0,e_bit),dec2bin(i+1,m)] ;
                        break;
                        
                    end
                end
            end
        end
        return;
        % if fraction >= (1/2+1/4+1/8+1/16)
        %     % 大于等于7/8到1的情况
        %     fp8Binary = [signBit '0001000'] ;
        % elseif fraction > 0.8125
        %     % 大于等于6/8到7/8的情况
        %     fp8Binary = [signBit '0000111'] ;
        % elseif fraction >= 0.6875
        %     % 大于等于5/8到6/8的情况
        %     fp8Binary = [signBit '0000110'] ;
        % elseif fraction > 0.5625
        %     % 大于等于4/8到5/8的情况
        %     fp8Binary = [signBit '0000101'] ;
        % elseif fraction >= 0.4375
        %     % 大于等于3/8到4/8的情况
        %     fp8Binary = [signBit '0000100'] ;
        % elseif fraction > 0.3125
        %     % 大于等于2/8到3/8的情况
        %     fp8Binary = [signBit '0000011'] ;
        % elseif fraction >= 0.1875   % 3/16
        %     % 大于等于1/8到2/8的情况
        %     fp8Binary = [signBit '0000010'] ;
        % elseif fraction > 0.0625  % 1/16
        %     % 大于等于0到1/8的情况
        %     fp8Binary = [signBit '0000001'] ;
        % else
        %     % 小于0.125的情况
        %     fp8Binary = [signBit '0000000'] ;
        % end
    else  % 规范数
        % 提取小数部分的前两位（四舍五入）
        fraction = fractionBits(1:m);
        threeBit = fractionBits(m);
        fourBit = fractionBits(m+1);
        % 判断 fractionBits 的第5位到最后一位是否存在1
        hasOnes = any(fractionBits(m+2:end) == '1');
        quan1 = repmat('1',1,m);
        if fourBit == '1'
            if isequal(fractionBits(1:m), quan1)
                fraction = repmat('0',1,m);
                exponent = exponent + 1;
            elseif hasOnes == 1 || threeBit == '1'
                fraction = dec2bin(bin2dec(fraction) + 1, m);
            end
        end
        % 构建fp8二进制表示
        fp8Binary = [signBit, dec2bin(exponent+bias, e_bit), fraction];
    end
end


