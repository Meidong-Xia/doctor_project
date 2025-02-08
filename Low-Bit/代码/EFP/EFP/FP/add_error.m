numSimulations = 1000;  % 每次累加的迭代次数
numIterations = 60;  % 累加次数

% 初始化误差存储
error8bit_e5m2 = zeros(numIterations, numSimulations);
error8bit_e4m3 = zeros(numIterations, numSimulations);
error8int = zeros(numIterations, numSimulations);
error16bit = zeros(numIterations, numSimulations);
error32bit = zeros(numIterations, numSimulations);

for sim = 1:numSimulations
    % 初始化累加变量
    sum64bit = 0;
    sum8bit_e5m2 = '00000000';
    sum8bit_e4m3 = '00000000';
    sum16bit = 0;
    sum32bit = 0;
    sum8int = 0;
    
    % 迭代累加
    for iter = 1:numIterations
        % 生成随机的8位整数
        % 生成小数点后10位的随机数
%         n = double(vpa(rand(), 10));
        n = 10*rand();
        
        % 8位加法 定点数
        value8int_bin = decimalToInt8(n);
        value8int_dec = int8ToDecimal(value8int_bin);
        sum8int = sum8int + value8int_dec;
        
        % 8位加法   e5m2
        value8bit_e5m2 = decimalTofp8_e5m2(n);
        sum8bit_e5m2 = fp8Add_e5m2(sum8bit_e5m2,value8bit_e5m2);
        sum8bit_dec_e5m2 = fp8Todecimal_e5m2(sum8bit_e5m2);
        
       % 8位加法    e4m3
        value8bit_e4m3 = decimalTofp8_e4m3(n);
        sum8bit_e4m3 = fp8Add_e4m3(sum8bit_e4m3,value8bit_e4m3);
        sum8bit_dec_e4m3 = fp8Todecimal_e4m3(sum8bit_e4m3);
        
        % 16位加法
        [~,~,sum16bit] = hp_add(sum16bit,n);

        % 32位加法
        sum32bit = single(single(sum32bit) + single(n));

        % 64位加法
        sum64bit = sum64bit + n;

        % 记录每次的累加误差
        error8bit_e5m2(iter, sim) = abs(double(sum64bit) - double(sum8bit_dec_e5m2))/abs(double(sum64bit));
        error8bit_e4m3(iter, sim) = abs(double(sum64bit) - double(sum8bit_dec_e4m3))/abs(double(sum64bit));
        error16bit(iter, sim) = abs(double(sum64bit) - double(sum16bit))/abs(double(sum64bit));
        error32bit(iter, sim) = abs(double(sum64bit) - sum32bit)/abs(double(sum64bit));
        error8int(iter, sim) = abs(double(sum64bit) - sum8int)/abs(double(sum64bit));
    end
end

% 计算每次累加的平均误差
averageError8bit_e5m2 = mean(error8bit_e5m2, 2);
averageError8bit_e4m3 = mean(error8bit_e4m3, 2);
averageError16bit = mean(error16bit, 2);
averageError32bit = mean(error32bit, 2);
averageError8int = mean(error8int, 2);
% 画图
figure;
semilogy(averageError8int, 'm', 'LineWidth', 2, 'DisplayName', 'INT 8');
hold on;
semilogy(averageError8bit_e5m2, 'r', 'LineWidth', 2, 'DisplayName', 'FP8 E5M2');
hold on;
semilogy(averageError8bit_e4m3, 'k', 'LineWidth', 2, 'DisplayName', 'FP8 E4M3');
hold on;
semilogy(averageError16bit, 'g', 'LineWidth', 2, 'DisplayName', 'FP16');
semilogy(averageError32bit, 'b', 'LineWidth', 2, 'DisplayName', 'FP32');
xlabel('累加次数');
ylabel('相对误差');
legend('show');
title('不同比特加法误差累积比较 随机数范围[0,10)');
grid on;
