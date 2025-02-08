% mex fp8Add_e5m2.cpp
% mex fp8Sub_e5m2.cpp
% mex i_hp_matrix_mul.cpp


%%  测试加法BASE 2 

    
    e_bit1 = 2;
    e_bit2 = 2;
    m_bit1 = 5; 
    m_bit2 = 5; 
    config1 = [e_bit1 m_bit1];
    config2 = [e_bit2 m_bit2];
    m = 2^(e_bit1+m_bit1+1);
    n = 2^(e_bit2+m_bit2+1);
    % Initialize result matrix
    mismatch_count = 0;
    mismatch_indices = [];
    differences = zeros(m * n, 1);
    data_pairs = cell(m * n, 7);
    
    
    tic;
    index = 1;
    
    for i = 1:m
        for j = 1:n
            % Generate binary strings
            fp8Binary1 = dec2bin(i - 1, e_bit1+m_bit1+1);
            fp8Binary2 = dec2bin(j - 1, e_bit2+m_bit2+1);
    
            % fp8Binary1 = '000000001110';
            % fp8Binary2 ='100000001111';
            % Call the MEX function
            [result_fp8_bin] = FP_add(fp8Binary1,fp8Binary2, config2);
             result_dec = fp8Todecimal(result_fp8_bin,config2);
    
    
            % Perform decimal addition for validation
            result_double_dec = fp8Todecimal(fp8Binary1, config1) + fp8Todecimal(fp8Binary2, config2) ;
            [result_double_fp8_bin]  = decimalTofp8(result_double_dec,config2);
             result_double_fp8_dec = fp8Todecimal(result_double_fp8_bin,config2);
            % Record data
            differences(index) = result_dec - result_double_fp8_dec;
            % differences(index) = bin2dec(result_fp8_bin)- bin2dec(result_double_fp8_bin);

            data_pairs{index, 1} = fp8Binary1;
            data_pairs{index, 2} = fp8Binary2;
            data_pairs{index, 3} = result_fp8_bin;
            data_pairs{index, 4} = result_double_fp8_bin;
            data_pairs{index, 5} = result_dec;
            data_pairs{index, 6} = result_double_fp8_dec;
            data_pairs{index, 7} = differences(index);
    
            % Compare results
            % if ~isequal(result_fp8_bin, result_double_fp8_bin)
            %     mismatch_count = mismatch_count + 1;
            %     mismatch_indices = [mismatch_indices; i, j];
            %             disp(fp8Binary1);
            % disp(fp8Binary2);
            % disp(result_fp8_bin);
            % % disp(result_bias);
            % disp(result_double_fp8_bin);
            % disp(result_dec);
            % disp(result_double_fp8_dec);
            % disp(differences(index));
            % break;
            % end
    
            if ~isequal(differences(index),0)&& ~isnan(differences(index))
            disp(fp8Binary1);
            disp(fp8Binary2);
            disp(result_fp8_bin);
            disp(result_double_fp8_bin);
            disp(result_dec);
            disp(result_double_fp8_dec);
            disp(differences(index));
            end
            % 
            index = index + 1;
        end
    end
    
    % Stop the timer
    elapsed_time = toc;
    
    % Display summary
    disp(['Total mismatches: ' num2str(mismatch_count)]);
    disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
    disp('Mismatched indices:');
    disp(mismatch_indices);
    
    % Plot differences
    figure;
    plot(1:(m * n), differences, '-o');
    xlabel('data');
    ylabel('Difference (fp8 - fp8->fp64->fp8)');
    title('Differences between fp8 mul and fp8->fp64->fp8 mul(e5m2)');
    
    % 假设 data_pairs 的第7列是要统计的列
    count_nonzero = 0;
    % 假设 data_pairs 的第7列是要统计的列
    filtered_data = {};
    
    % for i = 1:size(data_pairs, 1)
    %     % 检查第7位非零且不是NaN或Inf
    %     if data_pairs{i, 7} ~= 0 && ~isnan(data_pairs{i, 7}) && ~isinf(data_pairs{i, 7})
    %         count_nonzero = count_nonzero + 1;
    %                 % 将符合条件的行添加到新的数组中
    %         filtered_data = [filtered_data; data_pairs(i, :)];
    %     end
    % end
    
    disp(['第7位非零且不是NaN或Inf的个数: ' num2str(count_nonzero)]);


%%  测试减法

    
    e_bit1 = 2;
    e_bit2 = 2;
    m_bit1 = 5; 
    m_bit2 = 5; 
    config1 = [e_bit1 m_bit1];
    config2 = [e_bit2 m_bit2];
    m = 2^(e_bit1+m_bit1+1);
    n = 2^(e_bit2+m_bit2+1);
    % Initialize result matrix
    mismatch_count = 0;
    mismatch_indices = [];
    differences = zeros(m * n, 1);
    data_pairs = cell(m * n, 7);
    
    
    tic;
    index = 1;
    
    for i = 1:m
        for j = 1:n
            % Generate binary strings
            fp8Binary1 = dec2bin(i - 1, e_bit1+m_bit1+1);
            fp8Binary2 = dec2bin(j - 1, e_bit2+m_bit2+1);
    
            % fp8Binary1 = '000000001110';
            % fp8Binary2 ='100000001111';
            % Call the MEX function
            [result_fp8_bin] = FP_sub(fp8Binary1,fp8Binary2, config2);
             result_dec = fp8Todecimal(result_fp8_bin,config2);
    
    
            % Perform decimal addition for validation
            result_double_dec = fp8Todecimal(fp8Binary1, config1) - fp8Todecimal(fp8Binary2, config2) ;
            [result_double_fp8_bin]  = decimalTofp8(result_double_dec,config2);
             result_double_fp8_dec = fp8Todecimal(result_double_fp8_bin,config2);
            % Record data
            differences(index) = result_dec - result_double_fp8_dec;
            % differences(index) = bin2dec(result_fp8_bin)- bin2dec(result_double_fp8_bin);

            data_pairs{index, 1} = fp8Binary1;
            data_pairs{index, 2} = fp8Binary2;
            data_pairs{index, 3} = result_fp8_bin;
            data_pairs{index, 4} = result_double_fp8_bin;
            data_pairs{index, 5} = result_dec;
            data_pairs{index, 6} = result_double_fp8_dec;
            data_pairs{index, 7} = differences(index);
    
            % Compare results
            % if ~isequal(result_fp8_bin, result_double_fp8_bin)
            %     mismatch_count = mismatch_count + 1;
            %     mismatch_indices = [mismatch_indices; i, j];
            %             disp(fp8Binary1);
            % disp(fp8Binary2);
            % disp(result_fp8_bin);
            % % disp(result_bias);
            % disp(result_double_fp8_bin);
            % disp(result_dec);
            % disp(result_double_fp8_dec);
            % disp(differences(index));
            % break;
            % end
    
            if ~isequal(differences(index),0)&& ~isnan(differences(index))
            disp(fp8Binary1);
            disp(fp8Binary2);
            disp(result_fp8_bin);
            disp(result_double_fp8_bin);
            disp(result_dec);
            disp(result_double_fp8_dec);
            disp(differences(index));
            end
            % 
            index = index + 1;
        end
    end
    
    % Stop the timer
    elapsed_time = toc;
    
    % Display summary
    disp(['Total mismatches: ' num2str(mismatch_count)]);
    disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
    disp('Mismatched indices:');
    disp(mismatch_indices);
    
    % Plot differences
    figure;
    plot(1:(m * n), differences, '-o');
    xlabel('data');
    ylabel('Difference (fp8 - fp8->fp64->fp8)');
    title('Differences between fp8 mul and fp8->fp64->fp8 mul(e5m2)');
    
    % 假设 data_pairs 的第7列是要统计的列
    count_nonzero = 0;
    % 假设 data_pairs 的第7列是要统计的列
    filtered_data = {};
    
    % for i = 1:size(data_pairs, 1)
    %     % 检查第7位非零且不是NaN或Inf
    %     if data_pairs{i, 7} ~= 0 && ~isnan(data_pairs{i, 7}) && ~isinf(data_pairs{i, 7})
    %         count_nonzero = count_nonzero + 1;
    %                 % 将符合条件的行添加到新的数组中
    %         filtered_data = [filtered_data; data_pairs(i, :)];
    %     end
    % end
    
    disp(['第7位非零且不是NaN或Inf的个数: ' num2str(count_nonzero)]);


%%  e5m2加法
m = 256;
n = 256;

% Initialize result matrix
mismatch_count = 0;
mismatch_indices = [];
differences = zeros(m * n, 1);
data_pairs = cell(m * n, 7);

excel_filename = 'add_result.xlsx';

tic;
index = 1;

for i = 1:m
    for j = 1:n
        % Generate binary strings
        data1 = dec2bin(i - 1, 8);
        data2 = dec2bin(j - 1, 8);

        % Call the MEX function
        result_fp8_bin = fp8Add_e5m2(data1, data2);
        result_fp8_dec = fp8Todecimal_e5m2(result_fp8_bin);
        
        % Perform decimal addition for validation
        result_double_dec = fp8Todecimal_e5m2(data1) + fp8Todecimal_e5m2(data2);
        result_double_fp8_bin = decimalTofp8_e5m2(result_double_dec);
        result_double_fp8_dec = fp8Todecimal_e5m2(result_double_fp8_bin);
        % Record data
        differences(index) = result_fp8_dec - result_double_fp8_dec;
        
        data_pairs{index, 1} = data1;
        data_pairs{index, 2} = data2;
        data_pairs{index, 3} = result_fp8_bin;
        data_pairs{index, 4} = result_double_fp8_bin;
        data_pairs{index, 5} = result_fp8_dec;
        data_pairs{index, 6} = result_double_fp8_dec;
        data_pairs{index, 7} = differences(index);
        
        % Compare results
        if ~isequal(result_fp8_bin, result_double_fp8_bin)
            mismatch_count = mismatch_count + 1;
            mismatch_indices = [mismatch_indices; i, j];
                    disp(data1);
        disp(data2);
        disp(result_fp8_bin);
        disp(result_double_fp8_bin);
        disp(result_fp8_dec);
        disp(result_double_fp8_dec);
        disp(differences(index));
        end
        
        if ~isequal(differences(index),0)&& ~isnan(differences(index))
        disp(data1);
        disp(data2);
        disp(result_fp8_bin);
        disp(result_double_fp8_bin);
        disp(result_fp8_dec);
        disp(result_double_fp8_dec);
        disp(differences(index));
        end
        
        index = index + 1;
    end
end

% Stop the timer
elapsed_time = toc;

% Display summary
disp(['Total mismatches: ' num2str(mismatch_count)]);
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
disp('Mismatched indices:');
disp(mismatch_indices);

% Plot differences
figure;
plot(1:(m * n), differences, '-o');
xlabel('data');
ylabel('Difference (fp8 - fp8->fp64->fp8)');
title('Differences between fp8 Add and fp8->fp64->fp8 Add(e5m2)');

% 假设 data_pairs 的第7列是要统计的列
count_nonzero = 0;
% 假设 data_pairs 的第7列是要统计的列
filtered_data = {};

for i = 1:size(data_pairs, 1)
    % 检查第7位非零且不是NaN或Inf
    if data_pairs{i, 7} ~= 0 && ~isnan(data_pairs{i, 7}) && ~isinf(data_pairs{i, 7})
        count_nonzero = count_nonzero + 1;
                % 将符合条件的行添加到新的数组中
        filtered_data = [filtered_data; data_pairs(i, :)];
    end
end

disp(['第7位非零且不是NaN或Inf的个数: ' num2str(count_nonzero)]);

%%  e4m3加法
% mex fp8Add_e4m3.cpp
% mex fp8Sub_e4m3.cpp
m = 256;
n = 256;

% Initialize result matrix
mismatch_count = 0;
mismatch_indices = [];
differences = zeros(m * n, 1);
data_pairs = cell(m * n, 7);

excel_filename = 'add_result.xlsx';

tic;
index = 1;

for i = 1:m
    for j = 1:n
        % Generate binary strings
        data1 = dec2bin(i - 1, 8);
        data2 = dec2bin(j - 1, 8);

        % Call the MEX function
        result_fp8_bin = fp8Add_e4m3(data1, data2);
        result_fp8_dec = fp8Todecimal_e4m3(result_fp8_bin);
        
        % Perform decimal addition for validation
        result_double_dec = fp8Todecimal_e4m3(data1) + fp8Todecimal_e4m3(data2);
        result_double_fp8_bin = decimalTofp8_e4m3(result_double_dec);
        result_double_fp8_dec = fp8Todecimal_e4m3(result_double_fp8_bin);
        % Record data
        differences(index) = result_fp8_dec - result_double_fp8_dec;
        
        data_pairs{index, 1} = data1;
        data_pairs{index, 2} = data2;
        data_pairs{index, 3} = result_fp8_bin;
        data_pairs{index, 4} = result_double_fp8_bin;
        data_pairs{index, 5} = result_fp8_dec;
        data_pairs{index, 6} = result_double_fp8_dec;
        data_pairs{index, 7} = differences(index);
        
        % Compare results
        if ~isequal(result_fp8_bin, result_double_fp8_bin)
            mismatch_count = mismatch_count + 1;
            mismatch_indices = [mismatch_indices; i, j];
                    disp(data1);
        disp(data2);
        disp(result_fp8_bin);
        disp(result_double_fp8_bin);
        disp(result_fp8_dec);
        disp(result_double_fp8_dec);
        disp(differences(index));
        end
        
        if ~isequal(differences(index),0)&& ~isnan(differences(index))
        disp(data1);
        disp(data2);
        disp(result_fp8_bin);
        disp(result_double_fp8_bin);
        disp(result_fp8_dec);
        disp(result_double_fp8_dec);
        disp(differences(index));
        end
        
        index = index + 1;
    end
end

% Stop the timer
elapsed_time = toc;

% Display summary
disp(['Total mismatches: ' num2str(mismatch_count)]);
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
disp('Mismatched indices:');
disp(mismatch_indices);

% Plot differences
figure;
plot(1:(m * n), differences, '-o');
xlabel('data');
ylabel('Difference (fp8 - fp8->fp64->fp8)');
title('Differences between fp8 add and fp8->fp64->fp8 add(e4m3)');

% 假设 data_pairs 的第7列是要统计的列
count_nonzero = 0;
% 假设 data_pairs 的第7列是要统计的列
filtered_data = {};

for i = 1:size(data_pairs, 1)
    % 检查第7位非零且不是NaN或Inf
    if data_pairs{i, 7} ~= 0 && ~isnan(data_pairs{i, 7}) && ~isinf(data_pairs{i, 7})
        count_nonzero = count_nonzero + 1;
                % 将符合条件的行添加到新的数组中
        filtered_data = [filtered_data; data_pairs(i, :)];
    end
end


%%  
% e4m3减法
m = 256;
n = 256;

% Initialize result matrix
mismatch_count = 0;
mismatch_indices = [];
differences = zeros(m * n, 1);
data_pairs = cell(m * n, 7);

excel_filename = 'add_result.xlsx';

tic;
index = 1;

for i = 1:m
    for j = 1:n
        % Generate binary strings
        data1 = dec2bin(i - 1, 8);
        data2 = dec2bin(j - 1, 8);

        % Call the MEX function
        result_fp8_bin = fp8Sub_e4m3(data1, data2);
        result_fp8_dec = fp8Todecimal_e4m3(result_fp8_bin);
        
        % Perform decimal addition for validation
        result_double_dec = fp8Todecimal_e4m3(data1) - fp8Todecimal_e4m3(data2);
        result_double_fp8_bin = decimalTofp8_e4m3(result_double_dec);
        result_double_fp8_dec = fp8Todecimal_e4m3(result_double_fp8_bin);
        % Record data
        differences(index) = result_fp8_dec - result_double_fp8_dec;
        
        data_pairs{index, 1} = data1;
        data_pairs{index, 2} = data2;
        data_pairs{index, 3} = result_fp8_bin;
        data_pairs{index, 4} = result_double_fp8_bin;
        data_pairs{index, 5} = result_fp8_dec;
        data_pairs{index, 6} = result_double_fp8_dec;
        data_pairs{index, 7} = differences(index);
        
        % Compare results
        if ~isequal(result_fp8_bin, result_double_fp8_bin)
            mismatch_count = mismatch_count + 1;
            mismatch_indices = [mismatch_indices; i, j];
                    disp(data1);
        disp(data2);
        disp(result_fp8_bin);
        disp(result_double_fp8_bin);
        disp(result_fp8_dec);
        disp(result_double_fp8_dec);
        disp(differences(index));
        end
        
        if ~isequal(differences(index),0)&& ~isnan(differences(index))
        disp(data1);
        disp(data2);
        disp(result_fp8_bin);
        disp(result_double_fp8_bin);
        disp(result_fp8_dec);
        disp(result_double_fp8_dec);
        disp(differences(index));
        end
        
        index = index + 1;
    end
end

% Stop the timer
elapsed_time = toc;

% Display summary
disp(['Total mismatches: ' num2str(mismatch_count)]);
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
disp('Mismatched indices:');
disp(mismatch_indices);

% Plot differences
figure;
plot(1:(m * n), differences, '-o');
xlabel('data');
ylabel('Difference (fp8 - fp8->fp64->fp8)');
title('Differences between fp8 Sub and fp8->fp64->fp8 Sub(e4m3)');

% 假设 data_pairs 的第7列是要统计的列
count_nonzero = 0;
% 假设 data_pairs 的第7列是要统计的列
filtered_data = {};

for i = 1:size(data_pairs, 1)
    % 检查第7位非零且不是NaN或Inf
    if data_pairs{i, 7} ~= 0 && ~isnan(data_pairs{i, 7}) && ~isinf(data_pairs{i, 7})
        count_nonzero = count_nonzero + 1;
                % 将符合条件的行添加到新的数组中
        filtered_data = [filtered_data; data_pairs(i, :)];
    end
end

disp(['第7位非零且不是NaN或Inf的个数: ' num2str(count_nonzero)]);
disp(['第7位非零且不是NaN或Inf的个数: ' num2str(count_nonzero)]);