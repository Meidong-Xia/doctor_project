function [result_fp8Binary, result_config] = EFP_sqrt(fp8Binary1,config1)
    % % 检查输入的长度是否正确
            % fp8Binary1 ='01000000000';
    % %         fp8Binary2 = '111111111111110';
    % e_bit1 = 2;
    % % e_bit2 = 2;
    % m_bit1 = 9; 
    % % m_bit2 = 12; 
    % % result_m_bit = 12;
    % bias1 = 8;
    % % bias2 = 2;
    % base = 10;
    e_bit1 = config1(1);
    m_bit1 = config1(2);
    bias1 = config1(3);

    result_config = zeros(1,3);
    result_config(1) = e_bit1;
    result_config(2) = m_bit1;
    result_m_bit = m_bit1;

    if numel(fp8Binary1) ~= e_bit1 + m_bit1 + 1 
        error('输入的二进制字符串长度不正确。');
    end
    
    % 提取符号、指数和尾数
    signBit1 = str2double(fp8Binary1(1));
    if signBit1 == 1
        result_fp8Binary = ['1',dec2bin(0,result_m_bit+e_bit1)];
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    end
    exponentBits1 = bin2dec(fp8Binary1(2:e_bit1+1));
    % fractionBits1 = bin2dec(fp8Binary1(e_bit1+2:e_bit1+m_bit1+1));
    

    
    if strcmp(fp8Binary1,dec2bin(0,m_bit1+e_bit1+1)) 
        result_fp8Binary = dec2bin(0,1+e_bit1+result_m_bit);
        result_config(1) = e_bit1;
        result_config(3) = bias1;
        return;
    end

    exp = exponentBits1-bias1;
    bias = 2;
    fractionBits_new = ['0',fp8Binary1(e_bit1+2:e_bit1+m_bit1)];
    % fraction = bin2dec(fp8Binary1(e_bit1+2:e_bit1+m_bit1+1));
    if mod(exp,2)==0
        % if exp <= -4
            % if (fraction==0&&exponentBits1~=0) %|| fraction == 1
            %     fractionBits_new = [fp8Binary1(e_bit1+2:e_bit1+m_bit1),'1'];
            % end
        % end
        exp = exp/2;
        exp = exp + bias;
        
    else
        exp = (exp-1)/2;
        exp = exp + bias;
        fractionBits = bin2dec(fractionBits_new) + 2^(m_bit1-1);
        if fractionBits >= 2^(m_bit1)
            fractionBits_new = dec2bin(fractionBits - 2^(m_bit1),m_bit1);
            exp = exp+1;
        else
            fractionBits_new = dec2bin(fractionBits,m_bit1);
        end
    end

    if exp > 3
        bias = bias - (exp-3);
        exp = 3;
    elseif exp < 0
        bias = bias - exp;
        exp = 0;
    end

    result_config(3) = bias;
    result_config(1) = e_bit1;
    if exp == 0 && bin2dec(fractionBits_new) == 0
        exp = 1;
        result_config(3) = result_config(3) +1;
    end
    result_fp8Binary = [num2str(signBit1), dec2bin(exp, result_config(1)), fractionBits_new];
        
end
