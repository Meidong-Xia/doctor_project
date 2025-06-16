M = 64;
fc = 24*1e9;
F = 0.2*fc;
N = 1000;
S = 5;
gap1 = F/S;
gap2 = 2*pi/N;
bp = zeros(N,S+1);
v2 = f_calArraySteerVector(M,pi/3,0,fc);
j = 1;
for f = 0:gap1:F
    i = 1;
    for theta = gap2:gap2:2*pi
        v1 = f_calArraySteerVector(M,theta,f,fc);
        if theta > pi/2 
            bp(i,j) = -50;
        else
            bp(i,j) = 20*log10(abs(v1.'*conj(v2))^2);
            if bp(i,j) < -50
                bp(i,j) = -50;
            end
        end
        i = i+1;
    end
    j = j + 1;
end
color = ['yellow',]
for j = 1:S+1
    polarplot(gap2:gap2:2*pi,(bp(:,j)),'LineWidth', 1.5)
    hold on
end
% thetaticks([0  90 180 270 360]); % 设置刻度
thetalim([0,360]); % 限制角度范围

rlim([-50,0]); % 设置 dB 轴范围，例如 [-40, 0] dB
rticks([-50,-40,-30, -20, -10, 0]);
rticklabels([""," ","", "-20 dB", " ", "0 dB"])

ax = gca; % 获取当前坐标轴
ax.GridAlpha = 0.5; % 网格完全不透明
ax.LineWidth = 0.6; % 加粗极坐标的圆线

print(gcf, '波束斜视.png', '-dpng', '-r600'); % 600 DPI 高分辨率 PNG
