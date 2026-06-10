numTests = 5;
dimensions = [50 100 150 200 250 300];  % 从8维到256维，每次增加8维
sparsity = 0.5;  % 设置稀疏度为0.5
% methods = {'Gaussian inverse','cholesky inverse','CSM','improved Neumann'};
%  methods = {'Gaussian inverse'};
% colors_1 = {'b-s', 'r-s', 'g-s', 'm-s'};
% colors_2 = {'b-*', 'r-*', 'g-*', 'c-*'};
% colors_3 = {'b-o','r-o','g-o', 'r-o'};
% colors_4 = {'b-*', 'r-s', 'g-o', 'm-+'};
% numMethods = length(methods);
% average_errors_single = zeros(length(dimensions),1);
average_errors_double = zeros(length(dimensions),1);
average_errors_half = zeros(length(dimensions),1);
double mul_double;

% for methodIndex = 1:numMethods
    for dimIndex = 1:length(dimensions)
        dimension = dimensions(dimIndex);
        fprintf('当前维度：%d\n', dimension);
        error_sum_single = 0;
        error_sum_double = 0;
        error_sum_half = 0;
        
        for test = 1:numTests
            A = randn(dimension);
            B = randn(dimension);
            C = zeros(dimension);
%             mul_single = single(A)*single(B);
            mul_double = double(A)*double(B);
            [~, ~, mul_half] = r_hp_matrix_mul(A,B);
            [~, ~, mul_double] = r_hp_matrix_add(C,mul_double);
            mul=A*B;

            % 计算误差
%              error_single = norm(mul_single - mul, 'fro') ./ norm(mul, 'fro');
            error_double = norm(mul_double - mul, 'fro') ./ norm(mul, 'fro');
             error_half = norm(mul_half - mul, 'fro') ./ norm(mul, 'fro');

            % 累积误差
%              error_sum_single = error_sum_single + error_single;
            error_sum_double = error_sum_double + error_double;
             error_sum_half = error_sum_half + error_half;
        end

         % 计算平均误差
%          average_errors_single(dimIndex) = error_sum_single / numTests;
        average_errors_double(dimIndex) = error_sum_double / numTests;
         average_errors_half(dimIndex) = error_sum_half / numTests;
    end
% end

% semilogy(dimensions, average_errors_half, colors_4{1}, 'LineWidth', 1, 'MarkerSize', 6);
% hold on;
% semilogy(dimensions, average_errors_single, colors_4{2}, 'LineWidth', 1, 'MarkerSize', 6);
% semilogy(dimensions, average_errors_double, colors_4{3}, 'LineWidth', 1, 'MarkerSize', 6, 'HandleVisibility', 'on');
% 
% xlabel('维度');
% ylabel('相对误差 (对数坐标)');
% title('不同精度下矩阵乘法误差 (FP16、FP32、FP64)');
% grid on;
% 
% legend('FP16', 'FP32', 'FP64'); % 添加曲线标注
% hold off;




