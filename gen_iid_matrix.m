%% 生成与输入矩阵中的元素独立同分布的新三维随机矩阵
% 输入参数：
%   input_matrix - 输入矩阵，维数为1 2 3等均可，因为其统计特性和维数无关
%   rows - 新矩阵的行数
%   cols - 新矩阵的列数
%   pages - 新矩阵的页数
% 输出参数：
%   iid_matrix - 生成的独立同分布矩阵，维数为(rows, cols, pages)

function iid_matrix = gen_iid_matrix(input_matrix, rows, cols, pages)

    % 将输入矩阵展平为一维向量
    data = input_matrix(:);
    
    % 计算需要生成的总元素数
    total_elements = rows * cols * pages;
    
    % 从原始数据中随机抽样（有放回），保持独立同分布特性
    indices = randi(length(data), total_elements, 1);
    sampled_data = data(indices);
    
    % 重塑为目标维度
    iid_matrix = reshape(sampled_data, rows, cols, pages);
end
