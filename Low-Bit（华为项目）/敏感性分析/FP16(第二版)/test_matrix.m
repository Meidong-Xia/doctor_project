clc;
clear;
% 设置循环次数
numTests = 1000;
% 初始化误差累积变量
totalError = 0;
% 设置异常值检测阈值
threshold = 0.5;  % 你可以根据需要调整阈值
% 初始化异常值计数器
outlierCount = 0;
% 初始化历史记录数组
inputHistoryA = zeros(4, 4, numTests);
inputHistoryB = zeros(4, 4, numTests);
for i = 1:numTests
    % 生成随机4x4复数矩阵
%     a = rand(4, 4)+rand(4, 4) * 1i;
%     b = rand(4, 4) +rand(4, 4) * 1i;
    b = rand()+rand() * 1i;
% 
    a = rand(4, 4);
%     b = rand();
     a(1,1)=1;
    % 记录输入数据的历史
    inputHistoryA(:, :, i) = a;
    inputHistoryB(:, :, i) = b;
    % 调用MEX函数
    [a_halfprecision, b_halfprecision, sum_ab] = i_hp_matrix_mul_e(complex(a), complex(b));
    % 计算误差并累积
    true = a.*b;
    error = norm(sum_ab - (true))/norm(true);
    totalError = totalError + error;
    % 如果误差超出阈值，则记录为异常值
    if error > threshold
        outlierCount = outlierCount + 1;
        disp(['Outlier detected in test ', num2str(i), ', error: ', num2str(error)]);
        disp(['Input matrix a:']);
        disp(a);
        disp(['Input matrix b:']);
        disp(b);
        disp(['Result matrix sum_ab:']);
        disp(sum_ab);
        disp(['true:']);
        disp(true);
    end
end
% 计算平均误差
averageError = totalError / numTests;
% 显示结果
disp(['Average error over ', num2str(numTests), ' tests: ', num2str(averageError)]);
disp(['Total outliers: ', num2str(outlierCount)]);
