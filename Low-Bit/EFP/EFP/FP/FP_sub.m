function fp8Res = FP_sub(data1, data2,config)
    e_bit = config(1);
    m_bit = config(2);
    % e_bit = 2;
    % m_bit = 5;
    % % 判断是否为无效数：阶码全1，尾数非全0，则为无效数
    % data1 = '00000001';
    % data2 = '01011110';
    if numel(data1) ~= 1+e_bit+m_bit ||numel(data2) ~= 1+e_bit+m_bit
        error('输入的二进制字符串必须为8位。');
    end
        % 获取符号位、阶码和尾数
        sign1 = (data1(1));
        sign2 = dec2bin((1-bin2dec(data2(1))));
        exp1 = data1(2:e_bit+1);
        exp2 = data2(2:e_bit+1);
        mant1 = data1(e_bit+2:e_bit+m_bit+1);
        mant2 = data2(e_bit+2:e_bit+m_bit+1);

        % 1. 检查操作数中是否有0、Inf、NaN
        if  (bin2dec(exp1) == 2^e_bit-1 && bin2dec(mant1) == 2^m_bit-1) || (bin2dec(exp2) == 2^e_bit-1 && bin2dec(mant2) == 2^m_bit-1)
            fp8Res = ['0', repmat('1', 1, e_bit+m_bit)]; % NAN 和任何数相加为 NAN (0 1111 111)
            return;
        end

        if (bin2dec(exp1) == 0 && bin2dec(mant1) == 0)
            fp8Res = [sign2 data2(2:end)]; % 0 + x = x
            return;
        elseif (bin2dec(exp2) == 0 && bin2dec(mant2) == 0)
            fp8Res = data1;
            return;
        end

        mant1 = bin2dec(mant1);
        mant2 = bin2dec(mant2);
        % 判断是否为非规格化值，并进行规格化处理
        if bin2dec(exp1) == 0
            exp1 = '1'; % 非规格化值指数位为 1-b = 1-7 = -6
            mant1 = bitshift(mant1, 15); % 或0 代表没有前导数1
        else
            mant1 = bitshift(mant1+2^m_bit, 15); % 或0x08是加上前导数1，左移15位是为了对阶
        end

        if bin2dec(exp2) == 0
            exp2 = '1';
            mant2 = bitshift(mant2, 15);
        else
            mant2 = bitshift(mant2+2^m_bit, 15);
        end

        % 进行对阶操作，确定新exp，后续只需要进行尾数的累加
        exp1_dec = bin2dec(exp1);
        exp2_dec = bin2dec(exp2);
        if exp1_dec > exp2_dec
            expNew = exp1_dec;
            mant2 = bitshift(mant2, exp2_dec - exp1_dec);
            exp2_dec = exp1;
        else
            expNew = exp2_dec;
            mant1 = bitshift(mant1, exp1_dec - exp2_dec);
            exp1_dec = exp2_dec;
        end

        % 尾数累加
        if xor(bin2dec(sign1), bin2dec(sign2)) % 一正一负
            if mant1 > mant2
                mantNew = mant1 - mant2;
                signNew = sign1;
            else
                mantNew = mant2 - mant1;
                signNew = sign2;
            end
        else % 同正或同负
            mantNew = mant1 + mant2;
            signNew = sign1;
        end

        lastBit = 0;
        % 产生进位，并进行规格化
        mantToExp = bitshift(mantNew, -15-m_bit); % 取出尾数的最高位
        if bitand(mantToExp, hex2dec('02')) % 有进位
            expNew = expNew + 1;
            lastBit = bitand(mantNew, hex2dec('01'));
            mantNew = bitshift(mantNew, -1);
        elseif mantToExp == 0 % 没有前导数1
            if mantNew == 0 % 尾数全为0
                fp8Res = repmat('0', 1, 1+e_bit+m_bit);
                return;
            else
                for i = 1:15+m_bit % 规格化
                    expNew = expNew - 1;
                    mantNew = bitshift(mantNew, 1);
                    if bitand(bitshift(mantNew, -15-m_bit), hex2dec('01')) % 判断前导数是否为1
                        break;
                    end
                end
            end
        end

        % Round to nearest even
        threeBit = bitand(bitshift(mantNew, -15),hex2dec('01')); % 第3位尾数是否为1
        fourBit = bitand(bitshift(mantNew, -14),hex2dec('01')); % 第4位尾数是否为1
        lastBit = bitand(mantNew, (bitshift(1, 14) - 1)) | lastBit; % 后14是否存在1
        mantNew = bitshift(mantNew, -15);
        if fourBit && (threeBit || lastBit)
            mantNew = mantNew + 1; % 舍入
        end

        if bitand(mantNew, 2^(m_bit+1)) % 舍入之后产生进位
            expNew = expNew + 1;
            mantNew = bitshift(mantNew, -1);
        end

        % 溢出判断
        if (expNew == 2^e_bit - 1 && mantNew >= 2^(m_bit+1)-1) || (expNew > 2^e_bit - 1) % 上溢出
                fp8Res = [signNew,repmat('1', 1, e_bit+m_bit-1),'0']; % 返回448
                return;
        elseif expNew <= 0 
            mantNew = bitshift(mantNew, expNew-1);
            expNew = 0;
        elseif expNew <= -m_bit
            fp8Res = 0;
            return;
        end

        mant_new = dec2bin(mantNew,m_bit);
        if length(mant_new) > m_bit
            mant_new = mant_new(2:end);
        end
        % 构造最终结果
        fp8Res = [signNew,dec2bin(expNew,e_bit),mant_new];
end

