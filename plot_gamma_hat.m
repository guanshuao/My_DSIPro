function plot_gamma_hat(varargin)
% PLOT_gamma_HAT  绘制不同窗口大小下估计相关系数期望与真实相关系数的关系曲线
%
% 用法:
%   plot_gamma_hat(n1, n2, ...)
%
% 输入参数:
%   n1, n2, ... : 1到4个正整数，代表不同的窗口尺寸
%
% 功能:
%   利用 LUT_lookup('forward', ...) 查询指定窗口尺寸n和一系列gamma值对应的
%   估计相关系数均值和标准差，在多张图上绘制：
%   - 图1: E(gamma_hat) 关于 gamma 的曲线
%   - 图2: σ(gamma_hat) 关于 gamma 的曲线
%   - 图3: E(gamma_hat) ± σ(gamma_hat) 的不确定度带
%
% 示例:
%   plot_gamma_hat(5, 9, 15)

    % 检查输入参数个数
    if nargin < 1 || nargin > 4
        error('plot_gamma_hat:InputCount', '请输入1到4个整数作为窗口尺寸。');
    end

    % 验证输入是否为标量正整数
    window_sizes = zeros(1, nargin);
    for i = 1:nargin
        val = varargin{i};
        if ~isscalar(val) || ~isnumeric(val) || val <= 0 || val ~= floor(val)
            error('plot_gamma_hat:InvalidInput', '所有输入必须为正整数。');
        end
        window_sizes(i) = val;
    end

    % 设置字号大小参数参数
    fs_label = 19;   % 横纵轴名称的字号
    fs_legend = 19;  % 图例的字号
    fs_tick = 19;    % 横纵坐标数字的字号
    fs_title = 12;   % 标题的字号

    colors = {'r', 'g', 'b', 'y'};
    gamma_values = linspace(0, 0.997, 998);

    % 预分配数据结构
    data = struct('n', {}, 'mean', {}, 'std', {}, 'color', {});
    legend_str = cell(1, nargin);
    
    % 查询LUT数据
    for i = 1:nargin
        n = window_sizes(i);
        color = colors{i};
        
        try
            % 正向查找: 给定n和gamma，获取E(gamma_hat)和σ(gamma_hat)
            [m, s] = LUT_lookup('forward', n, gamma_values);
        catch ME
            error('plot_gamma_hat:LookupError', ...
                '调用LUT_lookup失败 (n=%d): %s', n, ME.message);
        end
        
        data(i).n = n;
        data(i).mean = m;
        data(i).std = s;
        data(i).color = color;
        legend_str{i} = sprintf('$n=%d$', n);
    end

    % 图1: E(gamma_hat) vs gamma
    figure(1); 
    clf; 
    set(gcf, 'Color', 'w');
    hold on; 
    grid on; 
    box on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 1);
    
    for i = 1:length(data)
        plot(gamma_values, data(i).mean, 'Color', data(i).color, 'LineWidth', 2);
    end
    
    plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', '$\gamma$');
    
    xlabel('$|\gamma|$', 'Interpreter', 'latex', 'FontSize', fs_label);
    ylabel('$\mathrm{E}(|\hat{\gamma}|)$', 'Interpreter', 'latex', 'FontSize', fs_label);
    % title('$\mathrm{E}(|\hat{\gamma}|)$ - $|\gamma|$ Relationship', 'Interpreter', 'latex', 'FontSize', fs_title);
    axis equal;
    axis square;
    xlim([-0.02 1.02]);
    ylim([-0.02 1.02]);
    legend([legend_str, {'Unbiased Estimator'}], 'Interpreter', 'latex', 'Location', 'Best', 'FontSize', fs_legend);
    hold off;

    % 图2: σ(gamma_hat) vs gamma
    figure(2); 
    clf; 
    set(gcf, 'Color', 'w');
    hold on; 
    grid on; 
    box on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 1);

    for i = 1:length(data)
        plot(gamma_values, data(i).std, 'Color', data(i).color, 'LineWidth', 2);
    end

    xlabel('$|\gamma|$', 'Interpreter', 'latex', 'FontSize', fs_label);
    ylabel('$\sigma_{|\hat{\gamma}|}$', 'Interpreter', 'latex', 'FontSize', fs_label);
    % title('$\sigma_{|\hat{\gamma}|}$ - $|\gamma|$ Relationship', 'Interpreter', 'latex', 'FontSize', fs_title);
    xlim([-0.02 1.05]);
    axis square;
    legend(legend_str, 'Interpreter', 'latex', 'Location', 'Best', 'FontSize', fs_legend);
    hold off;

    % 图3: E(gamma_hat) ± σ(gamma_hat) vs gamma
    figure(3); 
    clf; 
    set(gcf, 'Color', 'w');
    hold on; 
    grid on; 
    box on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 1);
    for i = 1:length(data)
        lower_bound = data(i).mean - data(i).std;
        upper_bound = data(i).mean + data(i).std;
        
        % 绘制不确定度带 (半透明，不加入图例)
        fill([gamma_values, fliplr(gamma_values)], ...
             [lower_bound, fliplr(upper_bound)], ...
             data(i).color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
             
        % 绘制均值曲线
        plot(gamma_values, data(i).mean, ...
            'Color', data(i).color, 'LineWidth', 2, 'DisplayName', legend_str{i});
    end

    % 绘制参考线 y=gamma
    plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'Unbiased Estimator');
    
    xlabel('$|\gamma|$', 'Interpreter', 'latex', 'FontSize', fs_label);
    ylabel('$\mathrm{E}(|\hat{\gamma}|)$', 'Interpreter', 'latex', 'FontSize', fs_label);
    % title('$\mathrm{E}(|\hat{\gamma}|) \pm \sigma_{|\hat{\gamma}|}$ - $|\gamma|$ Relationship', 'Interpreter', 'latex', 'FontSize', fs_title);
    axis equal;
    axis square;
    xlim([-0.02 1.02]);
    ylim([-0.02 1.02]);
    legend('Interpreter', 'latex', 'Location', 'Best', 'FontSize', fs_legend);
    hold off;

end
