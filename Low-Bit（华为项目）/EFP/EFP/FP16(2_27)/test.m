% mex -setup %如果没有选择编译器，则执行该命令 


% mex i_hp_sub.cpp %已测试ok     0.000606    %二轮已测试  % 第三轮测试:  0.00055  用m
% mex i_hp_add.cpp %已测试ok     0.000240    %二轮已测试   % 第三轮测试：0.000240  用c
% mex i_hp_mul.cpp %已测试ok     0.000361    %二轮已测试   % 第三轮测试:  0.000302  用m 
% mex i_hp_div.cpp %已测试ok     0.000453     %二轮已测试   % 第三轮测试: 0.000304  用m
% mex c_hp_sqrt.cpp %已测试ok       0.000447   %二轮已测试   % 第三轮测试: 0.000194   用m  0.53047 seconds
% mex r_hp_sqrt.cpp %已测试ok       0.00017791
% 集成到  i_hp_sqrt 实数 复数都可ok
% mex i_hp_matrix_add.cpp %已测试ok    0.000218  0.13782s    % 第三轮测试:0.000217  用m
% mex i_hp_matrix_sub.cpp %已测试ok    0.000427  0.14166s  % 第三轮测试:0.000427  用m
% mex i_hp_matrix_mul.cpp %已测试ok    0.000317  0.19475s  
% 第三轮测试:0.000213 0.17145  用m
% mex i_hp_matrix_mul_e.cpp %已测试ok  0.000346 0.17847s  
% 第三轮测试:0.000291  0.16594s  用m
% mex i_hp_matrix_div_e.cpp %已测试ok  0.000457  0.28649s    
% 第三轮测试:0.000292  0.20327s 用m
% mex hp_turn.cpp
% mex c_hp_matrix_turn.cpp

% clc;
% clear;
% 设置循环次数
% a = rand(4,4) + rand(4,4) * 1i;
% b = improve_Gaussian_inv_16(a);
% true = inv(a);
% error = norm(b - (true), 'fro')/norm(true, 'fro');
numTests = 10000;
% 初始化误差累积变量
totalError = 0;
% 设置异常值检测阈值
threshold = 0.001;  % 你可以根据需要调整阈值
% 初始化异常值计数器
outlierCount = 0;
% 初始化历史记录数组
inputHistoryA = zeros(numTests, 1);
inputHistoryB = zeros(numTests, 1);
% a = hp_turn(0.49*2^-24);

% 记录开始时间
tic;



for i = 1:numTests
    % 生成随机输入数据（如果需要）
    a = 10000*(rand() + rand() * 1i);
    b = rand() + rand() * 1i;
   
%      a = rand(); 
% 
%      b = rand();
    % 记录输入数据的历史
    inputHistoryA(i) = a;
    inputHistoryB(i) = b;
    % 调用MEX函数
%     [a_halfprecision, b_halfprecision, result] = i_hp_div(a, b);
%       [result] = i_hp_div(a, b);
%      [a_halfprecision, result] = i_hp_sqrt(a);
      [result] = i_hp_sqrt(a);
%     [a_halfprecision] = hp_turn(a);
%     true = a/b;
    true_1=sqrt(a);
    % 计算误差并累积
    error = norm(abs(result - true_1))/norm(true_1);
%      error = norm(a_halfprecision - true)/norm(true);
    totalError = totalError + error;
    % 如果误差超出阈值，则记录为异常值
    if error > threshold
        outlierCount = outlierCount + 1;
        disp(['Outlier detected in test ', num2str(i), ', error: ', num2str(error)]);
        disp(['Input a: ', num2str(a)]);
         disp(['Input b: ', num2str(b)]);
        disp(['result: ', num2str(result)]);
        disp(['true: ', num2str(true_1)]);
    end
end

% 记录结束时间
elapsedTime = toc;

% 计算平均误差
averageError = totalError / numTests;
% 显示结果
disp(['Average error over ', num2str(numTests), ' tests: ', num2str(averageError)]);
disp(['Total outliers: ', num2str(outlierCount)]);
disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);
