%% 加载
% 假设您已经在当前工作目录中存储了所有可能用到的分数表文件
% 设置全局变量 fraction_tables
addpath('EFP');
base = 2;
global fraction_tables;
base_values = [2, 10];  % 支持的基数
max_m_bit = 16;  % 最大的位数
fraction_tables = cell(length(base_values), max_m_bit);

% 加载并存储所有可能用到的分数表
for i = 1:length(base_values)
    base = base_values(i);
    for m_bit = 1:max_m_bit
        folder_name = 'fraction_tables';  % 存储分数表的文件夹
        filename = fullfile(folder_name, ['fraction_table_base_', num2str(base), '_m_bit_', num2str(m_bit), '.mat']);
        loaded_data = load(filename);
        fraction_tables{i, m_bit} = loaded_data.fraction_table;
    end
end

% % 示例：使用分数表
% selected_base = 2;
% selected_m_bit = 5;
% selected_fraction_table = fraction_tables{selected_base == base_values, selected_m_bit};

global M1B10_tables;
m = 1;
diff = 0;
M1B10_tables = cell(4, diff+1);
folder_name = 'M1B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M1B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M1B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M2B10_tables;
m = 2;
diff = 1;
M2B10_tables = cell(4, diff+1);
folder_name = 'M2B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M2B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M2B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M3B10_tables;
m = 3;
diff = 1;
M3B10_tables = cell(4, diff+1);
folder_name = 'M3B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M3B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M3B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M4B10_tables;
m = 4;
diff = 2;
M4B10_tables = cell(4, diff+1);
folder_name = 'M4B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M4B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M4B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end


global M5B10_tables;

M5B10_tables = cell(4, 3);
folder_name = 'M5B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:3
        for s = 1:2
        filename = fullfile(folder_name, ['M5B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M5B10_tables{s, e} = double(loaded_data.(sprintf('M5B10_%d_%d_%d', i-1, e-1,s-1)));
        else
            M5B10_tables{s+2, e} = double(loaded_data.(sprintf('M5B10_%d_%d_%d', i-1, e-1,s-1)));
        end
        end
    end
end

global M6B10_tables;

M6B10_tables = cell(4, 3);
folder_name = 'M6B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:3
        for s = 1:2
        filename = fullfile(folder_name, ['M6B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M6B10_tables{s, e} = double(loaded_data.(sprintf('M6B10_%d_%d_%d', i-1, e-1,s-1)));
        else
            M6B10_tables{s+2, e} = double(loaded_data.(sprintf('M6B10_%d_%d_%d', i-1, e-1,s-1)));
        end
        end
    end
end

global M7B10_tables;
m = 7;
diff = 3;
M7B10_tables = cell(4, diff+1);
folder_name = 'M7B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M7B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M7B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M8B10_tables;
m = 8;
diff = 3;
M8B10_tables = cell(4, diff+1);
folder_name = 'M8B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M8B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M8B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M9B10_tables;
m = 9;
diff = 3;
M9B10_tables = cell(4, diff+1);
folder_name = 'M9B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M9B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M9B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M10B10_tables;
m = 10;
diff = 3;
M10B10_tables = cell(4, diff+1);
folder_name = 'M10B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M10B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M10B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

global M11B10_tables;
m = 11;
diff = 4;
M11B10_tables = cell(4, diff+1);
folder_name = 'M11B10_table_s';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        for s = 1:2
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '_',num2str(s-1),'.mat']);
        loaded_data = load(filename);
        if i ==1
            M11B10_tables{s, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        else
            M11B10_tables{s+2, e} = double(loaded_data.(sprintf('M%dB10_%d_%d_%d',m, i-1, e-1,s-1)));
        end
        end
    end
end

% 定义一个单一的cell数组来存储所有表
global M_B10_tables;
M_B10_tables = cell(1, 12); % 初始化cell数组，大小为12

% 假设你的表是通过某种方式加载到 M1B10_table, M2B10_table, ..., M12B10_table 中

% 将所有表存储到cell数组中
M_B10_tables{1} = M1B10_tables;
M_B10_tables{2} = M2B10_tables;
M_B10_tables{3} = M3B10_tables;
M_B10_tables{4} = M4B10_tables;
M_B10_tables{5} = M5B10_tables;
M_B10_tables{6} = M6B10_tables;
M_B10_tables{7} = M7B10_tables;
M_B10_tables{8} = M8B10_tables;
M_B10_tables{9} = M9B10_tables;
M_B10_tables{10} = M10B10_tables;
M_B10_tables{11} = M11B10_tables;

global M1B10_table;
m = 1;
diff = 0;
M1B10_table = cell(2, diff+1);
folder_name = 'M1B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M1B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M2B10_table;
m = 2;
diff = 1;
M2B10_table = cell(2, diff+1);
folder_name = 'M2B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M2B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M3B10_table;
m = 3;
diff = 1;
M3B10_table = cell(2, diff+1);
folder_name = 'M3B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M3B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M4B10_table;

M4B10_table = cell(2, 3);
folder_name = 'M4B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:3
        filename = fullfile(folder_name, ['M4B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M4B10_table{i, e} = double(loaded_data.(sprintf('M4B10_%d_%d', i-1, e-1)));
    end
end


global M5B10_table;

M5B10_table = cell(2, 3);
folder_name = 'M5B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:3
        filename = fullfile(folder_name, ['M5B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M5B10_table{i, e} = double(loaded_data.(sprintf('M5B10_%d_%d', i-1, e-1)));
    end
end

global M6B10_table;

M6B10_table = cell(2, 3);
folder_name = 'M6B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:3
        filename = fullfile(folder_name, ['M6B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M6B10_table{i, e} = double(loaded_data.(sprintf('M6B10_%d_%d', i-1, e-1)));
    end
end

global M7B10_table;
m = 7;
diff = 3;
M7B10_table = cell(2, diff+1);
folder_name = 'M7B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M7B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M8B10_table;
m = 8;
diff = 3;
M8B10_table = cell(2, diff+1);
folder_name = 'M8B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M8B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M9B10_table;
m = 9;
diff = 3;
M9B10_table = cell(2, diff+1);
folder_name = 'M9B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M9B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M10B10_table;
m = 10;
diff = 3;
M10B10_table = cell(2, diff+1);
folder_name = 'M10B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M10B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end

global M11B10_table;
m = 11;
diff = 4;
M11B10_table = cell(2, diff+1);
folder_name = 'M11B10_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B10_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M11B10_table{i, e} = double(loaded_data.(sprintf('M%dB10_%d_%d', m,i-1, e-1)));
    end
end


% 定义一个单一的cell数组来存储所有表
global M_B10_table;
M_B10_table = cell(1, 12); % 初始化cell数组，大小为12

% 假设你的表是通过某种方式加载到 M1B10_table, M2B10_table, ..., M12B10_table 中

% 将所有表存储到cell数组中
M_B10_table{1} = M1B10_table;
M_B10_table{2} = M2B10_table;
M_B10_table{3} = M3B10_table;
M_B10_table{4} = M4B10_table;
M_B10_table{5} = M5B10_table;
M_B10_table{6} = M6B10_table;
M_B10_table{7} = M7B10_table;
M_B10_table{8} = M8B10_table;
M_B10_table{9} = M9B10_table;
M_B10_table{10} = M10B10_table;
M_B10_table{11} = M11B10_table;

%
global M1B2_table;
m = 1;
diff = 3;
base = 2;
M1B2_table = cell(2, diff+1);
folder_name = 'M1B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M1B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M2B2_table;
m = 2;
diff = 4;
M2B2_table = cell(2, diff+1);
folder_name = 'M2B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M2B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M3B2_table;
m = 3;
diff = 5;
M3B2_table = cell(2, diff+1);
folder_name = 'M3B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M3B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M4B2_table;
m = 4;
diff = 6;
M4B2_table = cell(2, diff+1);
folder_name = 'M4B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M4B2_table{i, e} =double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end


global M5B2_table;
m = 5;
diff = 7;
M5B2_table = cell(2, diff+1);
folder_name = 'M5B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M5B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M6B2_table;
m = 6;
diff = 8;
M6B2_table = cell(2, diff+1);
folder_name = 'M6B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M6B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M7B2_table;
m = 7;
diff = 9;
M7B2_table = cell(2, diff+1);
folder_name = 'M7B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M7B2_table{i, e} =double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M8B2_table;
m = 8;
diff = 10;
M8B2_table = cell(2, diff+1);
folder_name = 'M8B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M8B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M9B2_table;
m = 9;
diff = 11;
M9B2_table = cell(2, diff+1);
folder_name = 'M9B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M9B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M10B2_table;
m = 10;
diff = 12;
M10B2_table = cell(2, diff+1);
folder_name = 'M10B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M10B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end

global M11B2_table;
m = 11;
diff = 13;
M11B2_table = cell(2, diff+1);
folder_name = 'M11B2_table';  % 存储分数表的文件夹
% 加载并存储所有可能用到的分数表
for i = 1:2
    for e = 1:diff+1
        filename = fullfile(folder_name, ['M',num2str(m),'B',num2str(base),'_', num2str(i-1), '_', num2str(e-1),'_table', '.mat']);
        loaded_data = load(filename);
        M11B2_table{i, e} = double(loaded_data.(sprintf('M%dB%d_%d_%d', m,base,i-1, e-1)));
    end
end


% 定义一个单一的cell数组来存储所有表
global M_B2_table;
M_B2_table = cell(1, 12); % 初始化cell数组，大小为12

% 假设你的表是通过某种方式加载到 M1B10_table, M2B10_table, ..., M12B10_table 中

% 将所有表存储到cell数组中
M_B2_table{1} = M1B2_table;
M_B2_table{2} = M2B2_table;
M_B2_table{3} = M3B2_table;
M_B2_table{4} = M4B2_table;
M_B2_table{5} = M5B2_table;
M_B2_table{6} = M6B2_table;
M_B2_table{7} = M7B2_table;
M_B2_table{8} = M8B2_table;
M_B2_table{9} = M9B2_table;
M_B2_table{10} = M10B2_table;
M_B2_table{11} = M11B2_table;


table = [M_B2_table; M_B10_table];
fraction_tables_par = parallel.pool.Constant(fraction_tables);
table_par = parallel.pool.Constant(table);