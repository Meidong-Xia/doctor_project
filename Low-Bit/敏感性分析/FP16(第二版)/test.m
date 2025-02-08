% mex -setup %如果没有选择编译器，则执行该命令 
% mex r_hp_add.cpp     0.0002217
% mex r_hp_sub.cpp     0.000136
% mex r_hp_mul.cpp     6.8558e-05
% mex r_hp_div.cpp     0.00087084
% mex r_hp_sqrt.cpp     0.00011103
% mex r_hp_matrix_add.cpp     0.00090782
% mex r_hp_matrix_sub.cpp     0.00053843
% mex r_hp_matrix_mul.cpp     0.0012748
% mex r_hp_matrix_div_e.cpp   0.0052686
% mex r_hp_matrix_mul_e.cpp    0.00030877

% mex i_hp_sub.cpp %已测试     0.000213    %二轮已测试
% mex i_hp_add.cpp %已测试     0.000336    %二轮已测试
% mex i_hp_mul.cpp %已测试     0.000207    %二轮已测试
% mex i_hp_div.cpp %已测试     0.000596     %二轮已测试
% mex i_hp_sqrt.cpp %已测试       0.000359   %二轮已测试
% mex i_hp_matrix_add.cpp %已测试    0.001270 
% mex i_hp_matrix_sub.cpp %已测试    0.000748 
% mex i_hp_matrix_mul.cpp %已测试    0.0027     
% mex i_hp_matrix_mul_e.cpp %已测试  0.00078
% mex i_hp_matrix_div_e.cpp %已测试  0.0024   
% mex hp_turn.cpp
clc;
clear;
% 设置循环次数
numTests = 1000;
% 初始化误差累积变量
totalError = 0;
% 设置异常值检测阈值
threshold = 0.001;  % 你可以根据需要调整阈值
% 初始化异常值计数器
outlierCount = 0;
% 初始化历史记录数组
inputHistoryA = zeros(numTests, 1);
inputHistoryB = zeros(numTests, 1);
for i = 1:numTests
    % 生成随机输入数据（如果需要）
    a = rand() + rand() * 1i;
%     b = rand() + rand() * 1i;
%      a = rand();
     b = rand();
    % 记录输入数据的历史
    inputHistoryA(i) = a;
    inputHistoryB(i) = b;
    % 调用MEX函数
    [a_halfprecision, b_halfprecision, result] = i_hp_div(a, b);
%      [a_halfprecision, result] = i_hp_sqrt(a);
%     [a_halfprecision] = hp_turn(a);
    true = a/b;
%     true=sqrt(a);
    % 计算误差并累积
    error = norm(result - true)/norm(true);
%      error = norm(a_halfprecision - true)/norm(true);
    totalError = totalError + error;
    % 如果误差超出阈值，则记录为异常值
    if error > threshold
        outlierCount = outlierCount + 1;
        disp(['Outlier detected in test ', num2str(i), ', error: ', num2str(error)]);
        disp(['Input a: ', num2str(a)]);
         disp(['Input b: ', num2str(b)]);
        disp(['result: ', num2str(result)]);
        disp(['true: ', num2str(true)]);
    end
end
% 计算平均误差
averageError = totalError / numTests;
% 显示结果
disp(['Average error over ', num2str(numTests), ' tests: ', num2str(averageError)]);
disp(['Total outliers: ', num2str(outlierCount)]);

