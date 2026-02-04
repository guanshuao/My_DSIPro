function tip(h,varargin)
% TIP -- 将图形窗口平铺显示在屏幕上
%
%    重置所有图形窗口的位置和顺序，使它们覆盖整个屏幕。
%
%    TIP(H) 对句柄 H 中的图形窗口执行此操作。
%    TIP(H1,H2,...) 对指定的句柄执行此操作。
%    
%    Example:
%      TIP         对主窗口的所有图形进行排序。
%      TIP(10:100) 仅对句柄在 10 到 100 之间的图形窗口进行排序。
%    
%    See also TOP, TRAP
%    

%// $Revision: 1.4 $  $Date: 2001/09/28 14:24:33 $
%// Bert Kampes, 11-Dec-2000 insar toolbox, changed implementation

%%% 获取主窗口的子对象（即图形窗口），并进行排序
switch nargin
  case 0
    % 如果没有输入参数，获取根对象 (0) 的所有子对象（即所有打开的 figure）
    qch = sort(get(0,'children'));
  otherwise
    % 如果有输入参数，处理输入的句柄
    % lying 可能是一个自定义函数，用于将数组展平为行向量
    qch = lying(h(ishandle(h)));
    for q = 1:length(varargin)
      h   = varargin{q};
      h   = lying(h(ishandle(h)));
      qch = [qch,h];
    end
end
%qch = unique(qch);
%%%
qnumfigs = length(qch);
if (qnumfigs==0)
  warning('No figure windows found, exiting.');
  return;
elseif (qnumfigs>10)
  warning('Too many figure windows, exiting.');
  return;
end

%%% 准备重新调整显示布局
set(0,'DefaultFigureUnits','pixels')
qsc       = (get(0,'screensize')); % 获取屏幕尺寸
qscwidth  = qsc(3); % 屏幕宽度
qscheight = qsc(4); % 屏幕高度
qnY       = 1; % 默认行数
if (qnumfigs>=4) qnY = 2; end; % 如果图形数量大于等于4，则分为2行
qnX       = ceil(qnumfigs/qnY); % 计算列数
qsizeX    = floor(qscwidth/ qnX); % 每个图形窗口的宽度
qsizeY    = floor(qscheight/qnY); % 每个图形窗口的高度
%
qfig = 0;
% 循环设置每个图形窗口的位置
for qY = 1:qnY
  for qX = 1:qnX
    % 计算当前图形窗口的位置 [left bottom width height]
    % 注意：MATLAB 的坐标原点在左下角
    qpos    = [qsc(1)+(qX-1)*qsizeX, ...
	      qsc(2)+rem(qY*qsizeY,qscheight), ...
	      qsizeX, floor(0.85*qsizeY)]; % 高度稍微缩小一点 (0.85)，留出标题栏空间
    qfig    = qfig + 1;
    if (qfig>qnumfigs) break; end;
    qhandle = qch(qfig);
    if ishandle (qhandle)
      set(qhandle,'Position',qpos); % 设置位置
      figure(qhandle); % 将该窗口置顶显示
    end
  end
end