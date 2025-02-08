function [result_fp8Binary, result_config] = EFP_sub(fp8Binary1, config1, fp8Binary2, config2, base, fraction_tables,table)
    % % 检查输入的长度是否正确
            % fp8Binary1 = '11111111';
            % fp8Binary2 = '01010001';
    % e_bit1 = 2;
    % e_bit2 = 2;
    % m_bit1 = 12; 
    % m_bit2 = 12; 
    % result_m_bit = 12;
    % bias1 = 2;
    % bias2 = 2;
    % base = 10;
    % table = table.Value;
    e_bit1 = config1(1);
    m_bit1 = config1(2);
    bias1 = config1(3);
    e_bit2 = config2(1);
    m_bit2 = config2(2);

    result_config = zeros(1,3);
    result_config(1) = e_bit1;
    result_config(2) = m_bit1;
    result_m_bit = m_bit1;

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
    elseif  strcmp(fp8Binary2,dec2bin(0,m_bit2+e_bit2+1))
        % 如果其中一个操作数为 0，将另一个操作数作为结果
            result_fp8Binary = fp8Binary1(1:1+result_m_bit+e_bit1);
            result_config(1) = e_bit1;
            result_config(3) = bias1;
        if bin2dec(result_fp8Binary(2:3)) == 0 && bin2dec(result_fp8Binary(4:end)) == 0
            result_fp8Binary(2:3) = dec2bin(1,2);
            result_config(3) = result_config(3) + 1;
        end
        return;
    end

    fp8Binary2(1) = dec2bin(1-bin2dec(fp8Binary2(1)));
    [result_fp8Binary, result_config] = EFP_add(fp8Binary1, config1, fp8Binary2, config2, base, fraction_tables,table);
end
