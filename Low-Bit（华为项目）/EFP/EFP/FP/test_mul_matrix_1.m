clc;
% mex i_fp8Mul_e5m2.cpp
% mex fp8Mul_Matrix_e5m2.cpp
% mex i_hp_matrix_mul.cpp
%%  单个复数计算 e5m2
num = 1000;
% 初始化误差
errors = zeros(num, 1);

% 执行1000次复数乘法并计算误差
for i = 1:num
    % 生成随机复数
    b = rand() + rand() * 1i;
    c = rand() + rand() * 1i;
    
    % 复数乘法
    [fp8RealMatrix_1, fp8ImagMatrix_1] = complexToFP8Matrix_e5m2(b);
    [fp8RealMatrix_2, fp8ImagMatrix_2] = complexToFP8Matrix_e5m2(c);
    [fp8result_real, fp8result_imag] = i_fp8Mul_e5m2(fp8RealMatrix_1{1}, fp8ImagMatrix_1{1}, fp8RealMatrix_2{1}, fp8ImagMatrix_2{1});
    complexMatrix = fp8MatrixToComplex_e5m2({fp8result_real}, {fp8result_imag});
    result = b*c;
    
    % 计算误差
    error = norm(result - complexMatrix)/norm(result);
    errors(i) = error;
end

% 求平均误差
averageError = mean(errors);

disp(['平均误差: ' num2str(averageError)]);

%%  单个复数计算 e4m3
% mex i_fp8Mul_e4m3.cpp
num = 1000;
% 初始化误差
errors = zeros(num, 1);

% 执行1000次复数乘法并计算误差
for i = 1:num
    % 生成随机复数
    b = rand() + rand() * 1i;
    c = rand() + rand() * 1i;
    
    % 复数乘法
    [fp8RealMatrix_1, fp8ImagMatrix_1] = complexToFP8Matrix_e4m3(b);
    [fp8RealMatrix_2, fp8ImagMatrix_2] = complexToFP8Matrix_e4m3(c);
    [fp8result_real, fp8result_imag] = i_fp8Mul_e4m3(fp8RealMatrix_1{1}, fp8ImagMatrix_1{1}, fp8RealMatrix_2{1}, fp8ImagMatrix_2{1});
    complexMatrix = fp8MatrixToComplex_e4m3({fp8result_real}, {fp8result_imag});
    result = b*c;
    
    % 计算误差
    error = norm(result - complexMatrix)/norm(result);
    errors(i) = error;
end

% 求平均误差
averageError = mean(errors);

disp(['平均误差: ' num2str(averageError)]);
%%   复数矩阵计算  e5m2
clear;
num = 1000;
% 初始化误差
errors = zeros(num, 1);
% 执行1000次复数乘法并计算误差
for i = 1:num
    % 生成随机复数
b = rand(4,4) + rand(4,4) * 1i;
c = rand(4,4) + rand(4,4) * 1i;
    
    % 复数乘法
    result = fp8Mul_Matrix_e5m2(b,c);
    true = b*c;
    
    % 计算误差
    error = norm(result - true)/norm(true);
    errors(i) = error;
end
% 求平均误差
averageError = mean(errors);
disp(['平均误差: ' num2str(averageError)]);


%%   复数矩阵计算  e4m3
clear;
num = 1000;
% 初始化误差
errors = zeros(num, 1);
% 执行1000次复数乘法并计算误差
for i = 1:num
    % 生成随机复数
b = rand(4,4) + rand(4,4) * 1i;
c = rand(4,4) + rand(4,4) * 1i;
    
    % 复数乘法
    result = fp8Mul_Matrix_e4m3(b,c);
    true = b*c;
    
    % 计算误差
    error = norm(result - true)/norm(true);
    errors(i) = error;
end
% 求平均误差
averageError = mean(errors);
disp(['平均误差: ' num2str(averageError)]);
%%    复数矩阵误差比较
clc;

% mex i_fp8Mul_e5m2.cpp
% mex fp8Mul_Matrix_e5m2.cpp
% mex fp8Mul_Matrix_e4m3.cpp
% mex i_hp_matrix_mul.cpp

% 设定维度数组
dimensions = [4, 8, 16, 32];
% dimensions = [128];
num = 2;
methods = 4;
% 初始化误差矩阵
errors = zeros(length(dimensions), methods);

% 执行每个维度的矩阵乘法并计算误差
for d = 1:length(dimensions)
    % 当前维度
    dim = dimensions(d);
    
    % 初始化当前维度下的误差
    currentErrors = zeros(num, methods);

    % 执行1000次复数乘法并计算误差
    for i = 1:num
        % 生成随机复数
        b = rand(dim, dim) + rand(dim, dim) * 1i;
        c = rand(dim, dim) + rand(dim, dim) * 1i;
        result = b * c;
        
        % 显示当前维度和计算次数
        disp(['维度: ' num2str(dim) ', 计算次数: ' num2str(i)]);
        
        % fp8 乘法 e5m2
        fp8_result_e5m2 = fp8Mul_Matrix_e5m2(b,c);
        disp('fp8 e5m2完成' );
        % 计算误差
        currentErrors(i, 1) = norm(result - fp8_result_e5m2) / norm(result);
        
        % fp8 乘法 e4m3
        fp8_result_e4m3 = fp8Mul_Matrix_e4m3(b,c);
        disp('fp8 e4m3完成' );
        % 计算误差
        currentErrors(i, 2) = norm(result - fp8_result_e4m3) / norm(result);

        % 16bit复数乘法
        [~,~,result16bit] = i_hp_matrix_mul(b, c);
        % 计算误差
        disp('fp16完成' );
        currentErrors(i, 3) = norm(result16bit - result) / norm(result);

        % 32bit复数乘法
        result32bit = fp32Mul_Matrix(b, c);
%           result32bit = single(b)*single(c);
          disp('fp32完成' );
        % 计算误差
        currentErrors(i, 4) = norm(result32bit - result) / norm(result);
    end

    % 求平均误差
    averageErrors = mean(currentErrors);

    % 存储平均误差
    errors(d, :) = averageErrors;
%       errors(d, :) = currentErrors;
end

% 显示结果
disp('维度 | 8-bit e5m2相对误差 | 8-bit e4m3相对误差 | 16-bit 相对误差 | 32-bit 相对误差');
disp([dimensions', errors]);
% 画图
figure;
semilogy(dimensions, errors(:, 1), '-o', 'LineWidth', 2, 'DisplayName', '8-bit e5m2 Error');
hold on;
semilogy(dimensions, errors(:, 2), '-o', 'LineWidth', 2, 'DisplayName', '8-bit e4m3 Error');
hold on;
semilogy(dimensions, errors(:, 3), '-o', 'LineWidth', 2, 'DisplayName', '16-bit Error');
semilogy(dimensions, errors(:, 4), '-o', 'LineWidth', 2, 'DisplayName', '32-bit Error');
xlabel('维度');
ylabel('相对误差');
legend('show');
title('不同比特位矩阵乘法误差比较');
grid on;
% 将 x 轴刻度设置为对数刻度
set(gca, 'YScale', 'log');
% 设置 x 轴刻度标签
xticks(dimensions);
xticklabels(cellstr(num2str(dimensions')));

%%
clc;
tic;
% mex i_fp8Mul_e5m2.cpp
% mex fp8Mul_Matrix_e5m2.cpp
% mex fp8Mul_Matrix_e4m3.cpp
% mex i_hp_matrix_mul.cpp

% 设定维度数组
dimensions = [4, 8, 16, 32, 64,128,256];
num = 2;
methods = 4;

% 初始化误差矩阵
errors = zeros(length(dimensions), methods);

% 执行每个维度的矩阵乘法并计算误差
for d = 1:length(dimensions)
    % 当前维度
    dim = dimensions(d);
    
    % 初始化当前维度下的误差
    currentErrors = zeros(num, methods);

    % 执行 num 次复数乘法并计算误差
    for i = 1:num
        % 生成随机复数
        b = rand(dim, dim) + rand(dim, dim) * 1i;
        c = rand(dim, dim) + rand(dim, dim) * 1i;
        result = b * c;
        
        % 显示当前维度和计算次数
        disp(['维度: ' num2str(dim) ', 计算次数: ' num2str(i)]);
        
        % fp8 乘法 e5m2
        fp8_result_e5m2 = fp8Mul_Matrix_e5m2(b,c);
        disp('fp8 e5m2完成' );
        % 计算误差
        currentErrors(i, 1) = norm(result - fp8_result_e5m2) / norm(result);
        
        % fp8 乘法 e4m3
        fp8_result_e4m3 = fp8Mul_Matrix_e4m3(b,c);
        disp('fp8 e4m3完成' );
        % 计算误差
        currentErrors(i, 2) = norm(result - fp8_result_e4m3) / norm(result);

        % 16bit复数乘法
        [~,~,result16bit] = i_hp_matrix_mul(b, c);
        disp('fp16 完成' );
        % 计算误差
        currentErrors(i, 3) = norm(result16bit - result) / norm(result);

        % 32bit复数乘法
        result32bit = fp32Mul_Matrix(b, c);
        disp('fp32 完成' );
        % 计算误差
        currentErrors(i, 4) = norm(result32bit - result) / norm(result);
    end

    % 求平均误差
    averageErrors = mean(currentErrors);

    % 存储平均误差
    errors(d, :) = averageErrors;
    
    % 存储每一次的误差
    allErrors{d} = currentErrors;
end

        % 结束计时并显示时间
        elapsedTime = toc;
        disp(['计算时间: ' num2str(elapsedTime) ' 秒']);
% 显示结果
disp('维度 | 8-bit e5m2相对误差 | 8-bit e4m3相对误差 | 16-bit 相对误差 | 32-bit 相对误差');
disp([dimensions', errors]);

% 画图
figure;
semilogy(dimensions, errors(:, 1), '-o', 'LineWidth', 2, 'DisplayName', 'FP8 E5M2');
hold on;
semilogy(dimensions, errors(:, 2), '-o', 'LineWidth', 2, 'DisplayName', 'FP8 E4M3');
hold on;
semilogy(dimensions, errors(:, 3), '-o', 'LineWidth', 2, 'DisplayName', 'FP16');
semilogy(dimensions, errors(:, 4), '-o', 'LineWidth', 2, 'DisplayName', 'FP32');
xlabel('维度');
ylabel('相对误差');
legend('show');
title('不同比特位矩阵乘法误差比较');
grid on;
% 将 x 轴刻度设置为对数刻度
set(gca, 'YScale', 'log');
% 设置 x 轴刻度标签
xticks(dimensions);
xticklabels(cellstr(num2str(dimensions')));
% 保存每一次的误差
save('allErrors_1.mat', 'allErrors');
