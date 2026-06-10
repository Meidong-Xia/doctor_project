% clc;
% clear;
% 设置循环次数
numTests = 1000;
% 初始化误差累积变量
totalError = 0;
% 设置异常值检测阈值
threshold = 0.5;  % 你可以根据需要调整阈值
% 初始化异常值计数器
outlierCount = 0;
% 初始化历史记录数组
di = 4;
inputHistoryA = zeros(di,di, numTests);
inputHistoryB = zeros(di,di, numTests);

% 记录开始时间
tic;

for i = 1:numTests
    % 生成随机4x4复数矩阵
    a = rand(di,di)+rand(di,di) * 1i;
%     b = rand(di,di) +rand(di,di) * 1i;
    b = rand()+rand() * 1i;
% 
%     a = rand(di, di);
     b = rand();
%      a(1,1)=1;
    % 记录输入数据的历史
    inputHistoryA(:, :, i) = a;
    inputHistoryB(:, :, i) = b;
    % 调用MEX函数
%     [a_halfprecision, b_halfprecision, sum_ab] = i_hp_matrix_div_e(complex(a), complex(b));
     sum_ab = i_hp_matrix_div_e(complex(double(a)), complex(double(b)));
    % 计算误差并累积
    true = a./b;
    error = norm(sum_ab - (true), 'fro')/norm(true, 'fro');
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

% 记录结束时间
elapsedTime = toc;

% 计算平均误差
averageError = totalError / numTests;
% 显示结果
disp(['Average error over ', num2str(numTests), ' tests: ', num2str(averageError)]);
disp(['Total outliers: ', num2str(outlierCount)]);
disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);