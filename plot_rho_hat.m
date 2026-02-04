function plot_rho_hat(varargin)
% PLOT_RHO_HAT  绘制不同窗口大小下估计相关系数期望与真实相关系数的关系曲线
%
% 用法:
%   plot_rho_hat(n1, n2, ...)
%
% 输入参数:
%   n1, n2, ... : 1到4个正整数，代表不同的窗口尺寸
%
% 功能:
%   利用 LUT_lookup('forward', ...) 查询指定窗口尺寸n和一系列rho值对应的
%   估计相关系数均值和标准差，在多张图上绘制：
%   - 图1: E(rho_hat) 关于 rho 的曲线
%   - 图2: σ(rho_hat) 关于 rho 的曲线
%   - 图3: E(rho_hat) ± σ(rho_hat) 的不确定度带
%
% 示例:
%   plot_rho_hat(5, 9, 15)

    % 检查输入参数个数
    if nargin < 1 || nargin > 4
        error('plot_rho_hat:InputCount', '请输入1到4个整数作为窗口尺寸。');
    end

    % 验证输入是否为标量正整数
    window_sizes = zeros(1, nargin);
    for i = 1:nargin
        val = varargin{i};
        if ~isscalar(val) || ~isnumeric(val) || val <= 0 || val ~= floor(val)
            error('plot_rho_hat:InvalidInput', '所有输入必须为正整数。');
        end
        window_sizes(i) = val;
    end

    colors = {'r', 'g', 'b', 'y'};
    rho_values = linspace(0, 0.997, 998);

    % 预分配数据结构
    data = struct('n', {}, 'mean', {}, 'std', {}, 'color', {});
    legend_str = cell(1, nargin);
    
    % 查询LUT数据
    for i = 1:nargin
        n = window_sizes(i);
        color = colors{i};
        
        try
            % 正向查找: 给定n和rho，获取E(rho_hat)和σ(rho_hat)
            [m, s] = LUT_lookup('forward', n, rho_values);
        catch ME
            error('plot_rho_hat:LookupError', ...
                '调用LUT_lookup失败 (n=%d): %s', n, ME.message);
        end
        
        data(i).n = n;
        data(i).mean = m;
        data(i).std = s;
        data(i).color = color;
        legend_str{i} = sprintf('n=%d', n);
    end

    % 图1: E(rho_hat) vs rho
    figure(1); 
    clf; 
    hold on; 
    grid on; 
    box on;
    
    for i = 1:length(data)
        plot(rho_values, data(i).mean, 'Color', data(i).color, 'LineWidth', 2);
    end
    
    plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', '\rho');
    
    xlabel('\rho', 'FontSize', 12);
    ylabel('E(\hat{\rho})', 'FontSize', 12);
    title('E(\hat{\rho}) - \rho Relationship', 'FontSize', 13);
    axis equal;
    axis square;
    xlim([0 1]);
    ylim([0 1]);
    legend([legend_str, {'\rho'}], 'Location', 'Best', 'FontSize', 10);
    hold off;

    % 图2: σ(rho_hat) vs rho
    figure(2); 
    clf; 
    hold on; 
    grid on; 
    box on;

    for i = 1:length(data)
        plot(rho_values, data(i).std, 'Color', data(i).color, 'LineWidth', 2);
    end

    xlabel('\rho', 'FontSize', 12);
    ylabel('\sigma_{\hat{\rho}}', 'FontSize', 12);
    title('\sigma_{\hat{\rho}} - \rho Relationship', 'FontSize', 13);
    xlim([0 1]);
    legend(legend_str, 'Location', 'Best', 'FontSize', 10);
    hold off;

    % 图3: E(rho_hat) ± σ(rho_hat) vs rho
    figure(3); 
    clf; 
    hold on; 
    grid on; 
    box on;
    for i = 1:length(data)
        lower_bound = data(i).mean - data(i).std;
        upper_bound = data(i).mean + data(i).std;
        
        % 绘制不确定度带 (半透明，不加入图例)
        fill([rho_values, fliplr(rho_values)], ...
             [lower_bound, fliplr(upper_bound)], ...
             data(i).color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
             
        % 绘制均值曲线
        plot(rho_values, data(i).mean, ...
            'Color', data(i).color, 'LineWidth', 2, 'DisplayName', legend_str{i});
    end

    % 绘制参考线 y=rho
    plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', '\rho');
    
    xlabel('\rho', 'FontSize', 12);
    ylabel('E(\hat{\rho})', 'FontSize', 12);
    title('E(\hat{\rho}) \pm \sigma_{\hat{\rho}} - \rho Relationship', 'FontSize', 13);
    axis equal;
    axis square;
    xlim([0 1]);
    ylim([0 1]);
    legend('Location', 'Best', 'FontSize', 10);
    hold off;

end
