function [result_fp8Binary, result_config] = EFP_add(fp8Binary1, config1, fp8Binary2, config2,base, fraction_tables,table)
    % % 检查输入的长度是否正确
    % fp8Binary1 ='00000001';
    % fp8Binary2 ='10000010';
    % e_bit1 = 2;
    % e_bit2 = 2;
    % m_bit1 = 5; 
    % m_bit2 = 5; 
    % bias1 = 2;
    % bias2 = 2;
    % config1 = [e_bit1 m_bit1 bias1];
    % config2 = [e_bit2 m_bit2 bias2];
    table = table.Value;

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
        error('输入的二进制字符串长度不正确。');
    end

    % 判断是否为NaN或者乘0
    if strcmp(fp8Binary1,['1',dec2bin(0,m_bit1+e_bit1)]) || strcmp(fp8Binary2,['1',dec2bin(0,m_bit2+e_bit2)])
        result_fp8Binary = ['1',dec2bin(0,result_m_bit+e_bit1)];
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    elseif strcmp(fp8Binary1,dec2bin(0,m_bit1+e_bit1+1)) && strcmp(fp8Binary2,dec2bin(0,m_bit2+e_bit2+1))
        result_fp8Binary = ['0',dec2bin(0,result_m_bit+e_bit1)];
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    elseif strcmp(fp8Binary1,dec2bin(0,m_bit1+e_bit1+1)) || strcmp(fp8Binary2,dec2bin(0,m_bit2+e_bit2+1))
        % 如果其中一个操作数为 0，将另一个操作数作为结果
        if strcmp(fp8Binary1, dec2bin(0, m_bit1 + e_bit1 + 1))
            result_fp8Binary = fp8Binary2(1:1+result_m_bit+e_bit2);
            result_config(1) = e_bit2;
            result_config(3) = bias2;
        else
            result_fp8Binary = fp8Binary1(1:1+result_m_bit+e_bit1);
            result_config(1) = e_bit1;
            result_config(3) = bias1;
        end
        if bin2dec(result_fp8Binary(2:3)) == 0 && bin2dec(result_fp8Binary(3:end)) == 0
            result_fp8Binary(2:3) = dec2bin(1,2);
            result_config(3) = result_config(3) +1;
        end
        return;
    end

    % 提取符号、指数和尾数
    signBit1 = str2double(fp8Binary1(1));
    exponentBits1 = bin2dec(fp8Binary1(2:e_bit1+1));
    fractionBits1 = bin2dec(fp8Binary1(e_bit1+2:e_bit1+m_bit1+1));
    
    signBit2 = str2double(fp8Binary2(1));
    exponentBits2 = bin2dec(fp8Binary2(2:e_bit2+1));
    fractionBits2 = bin2dec(fp8Binary2(e_bit2+2:e_bit2+m_bit2+1));

    % 计算指数和尾数的值
    exponent1 = exponentBits1 - bias1;
    exponent2 = exponentBits2 - bias2;
    e_diff = exponent1 - exponent2;

    max_m_bit = max(m_bit1,m_bit2);
    if max_m_bit < result_m_bit
        fractionBits1 = fractionBits1 * 2^(result_m_bit - m_bit1);
        fractionBits2 = fractionBits2 * 2^(result_m_bit - m_bit2);
        table_m_bit = result_m_bit;
    else
        fractionBits1 = fractionBits1 * 2^(max_m_bit - m_bit1);
        fractionBits2 = fractionBits2 * 2^(max_m_bit - m_bit2);
        table_m_bit = max_m_bit;
    end
    
    if e_diff < 0
        temp_fp8Binary = fp8Binary1;
        temp_e_bit = e_bit1;
        temp_m_bit = m_bit1;
        temp_bias = bias1;

        fp8Binary1 = fp8Binary2;
        e_bit1 = e_bit2;
        m_bit1 = m_bit2;
        bias1 = bias2;

        fp8Binary2 = temp_fp8Binary;
        e_bit2 = temp_e_bit;
        m_bit2 = temp_m_bit;
        bias2 = temp_bias;

        % 提取符号、指数和尾数
        signBit1 = str2double(fp8Binary1(1));
        exponentBits1 = bin2dec(fp8Binary1(2:e_bit1+1));
        fractionBits1 = bin2dec(fp8Binary1(e_bit1+2:e_bit1+m_bit1+1));
    
        signBit2 = str2double(fp8Binary2(1));
        exponentBits2 = bin2dec(fp8Binary2(2:e_bit2+1));
        fractionBits2 = bin2dec(fp8Binary2(e_bit2+2:e_bit2+m_bit2+1));
        % 计算指数和尾数的值
        exponent1 = exponentBits1 - bias1;
        exponent2 = exponentBits2 - bias2;
        e_diff = exponent1 - exponent2;
    end

            % 比较两个数的符号位
        if signBit1 == signBit2
            sign_result = signBit1;
        else
            % 如果符号位不同，则比较指数位
            if exponent1 == exponent2
                % 如果指数位相同，则比较小数位
                if fractionBits1 == fractionBits2
                    sign_result = 0;
                elseif fractionBits1 > fractionBits2
                    sign_result = signBit1;
                else
                    sign_result = signBit2;
                end
            elseif exponent1 > exponent2
                sign_result = signBit1;
            else
                sign_result = signBit2;
            end
        end
    % global base;
    if base == 10
        diff = [0 1 1 2 2 2 3 3 3 3 4 4];
    else
        diff = [3 4 5 6 7 8 9 10 11 12 13 14];
    end

    % 指数位之差大于2
    if abs(e_diff) > diff(table_m_bit)
        result_fp8Binary = fp8Binary1(1:1+result_m_bit+e_bit1);
            result_config(1) = e_bit1;
            result_config(3) = bias1;
        if bin2dec(result_fp8Binary(2:3)) == 0 && bin2dec(result_fp8Binary(3:end)) == 0
              result_fp8Binary(2:3) = dec2bin(1,2);
              result_config(3) = result_config(3) +1;
        end
        return;
    end

    %指数位之差小于等于2 确定查找表
    sign = xor(signBit1,signBit2);
    % 定义一个单一的cell数组来存储所有表

    % global M_B10_table;
    % global M_B2_table;
    if base == 10
        M_B10_table = table(2,:);
        table_all = M_B10_table{table_m_bit};
    else
        M_B2_table = table(1,:);
        table_all = M_B2_table{table_m_bit};
    end

    % if table_m_bit==1
    %     global M1B10_table;
    %     table_all = M1B10_table;
    % elseif table_m_bit==2
    %     global M2B10_table;
    %     table_all = M2B10_table;
    % elseif table_m_bit==3
    %     global M3B10_table;
    %     table_all = M3B10_table;
    % elseif table_m_bit==4
    %     global M4B10_table;
    %     table_all = M4B10_table;
    % elseif table_m_bit==5
    %     global M5B10_table;
    %     table_all = M5B10_table;
    % elseif table_m_bit==6
    %     global M6B10_table;
    %     table_all = M6B10_table;
    % elseif table_m_bit==7
    %     global M7B10_table;
    %     table_all = M7B10_table;
    % elseif table_m_bit==8
    %     global M8B10_table;
    %     table_all = M8B10_table;
    % elseif table_m_bit==9
    %     global M9B10_table;
    %     table_all = M9B10_table;
    % elseif table_m_bit==10
    %     global M10B10_table;
    %     table_all = M10B10_table;
    % elseif table_m_bit==11
    %     global M11B10_table;
    % global M_B10_table;
    %     table_all = M11B2_table;

    % elseif table_m_bit==12
    %     global M12B10_table;
    %     table_all = M12B10_table;
    % end


    % table_all = [M1B10_tables M2B10_tables M3B10_tables M4B10_tables M5B10_tables M6B10_tables M7B10_tables];

    % table_diff = 0;
    % for i = 1:table_m_bit-1
    %         table_diff = table_diff + diff(i)+1;
    % end
    % table = table_all{sign+1,e_diff+1+table_diff};
    table = table_all{sign+1,e_diff+1};

    if sign == 1 && e_diff == 0 && fractionBits1 == fractionBits2
        result_fp8Binary = dec2bin(0,result_m_bit+e_bit1+1);
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    end

    result_config(3) = bias1;

        if base == 10
            e_table_bit = 2;
            e_bias = 3;
        else
            e_table_bit = 4;
            e_bias = 15;
        end

        if fractionBits1 <= fractionBits2
            result_bin = dec2bin(table(1,fractionBits2-fractionBits1+1)+fractionBits1,e_table_bit+max_m_bit);
        else
            result_bin = dec2bin(table(2,fractionBits1-fractionBits2+1)+fractionBits2,e_table_bit+max_m_bit);
        end
        
             exponent_binary = bin2dec(result_bin(1:e_table_bit));
             fraction_binary = bin2dec(result_bin(e_table_bit+1:end));


            if sign == 0
                exponent_sum = exponent1 + exponent_binary;
                exponent_bit = exponent_sum + bias1;
                if exponent_bit > 3
                    result_config(3) = bias1 - (exponent_bit-3);
                    exponent_bit = 3;
                end
            else
                exponent_sum = exponent1 + exponent_binary-e_bias;
                exponent_bit = exponent_sum + bias1;
                if exponent_bit > 3
                    result_config(3) = bias1 - (exponent_bit-3);
                    exponent_bit = 3;
                elseif exponent_bit < 0 
                    result_config(3) = bias1 - (exponent_bit-0);
                    exponent_bit = 0;
                    if fraction_binary == 0
                        result_config(3) = result_config(3) + 1;
                        exponent_bit = 1;
                    end
                elseif exponent_bit == 0
                     if fraction_binary == 0
                        result_config(3) = result_config(3) + 1;
                        exponent_bit = 1;
                    end
                end
            end

        if exponent_bit == 0 && fraction_binary == 0
            exponent_bit = 1;
            result_config(3) = result_config(3) +1;
        end

        result_config(1) = e_bit1;
            % 返回乘积结果
        result_fp8Binary = [num2str(sign_result), dec2bin(exponent_bit, result_config(1)), dec2bin(fraction_binary, result_m_bit)];
        
end
