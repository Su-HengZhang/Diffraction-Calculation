close all;
clear,clc;
%%
%----load 物体采样矩阵----%
load('usaf4um.mat');
forg=usaf4um;
[R,~]=size(forg);

%----实际长度以mm为单位-----%
lambda=5e-4;    % 光波长
z1=80;          % 物体到透镜的距离
id=4e-3;        % 理想采样间隔

%%
M=round((lambda*z1)/id^2); % 第1次衍射计算时，计算窗口内的采样点数

u0=pad2center(forg,M,M); % 对物体进行零填充至 M*M

% 计算从物面到透镜前表面的衍射
u1=fresnelas(u0); % 利用角谱法计算
% u1=fresnelsft(u0); % 利用单次傅里叶变换算法计算

%%
lens_F=40;       % 透镜焦距
Q=round((lambda*lens_F)/id^2); % 透镜的无量纲化焦距 P*P 的矩形孔径

lens_dia=lens_F;    % 透镜的孔径为边长是lens_dia的矩形，需要 lens_dia<=lens_F
% 否则会导致对于透镜的二次相因子欠采样
P=round((lambda*lens_dia)/id^2);

x=0:(P-1);
x=x-floor(P/2);
y=x;
[Y,X]=meshgrid(y,x);

lens=exp(-1i*pi*(X.^2+Y.^2)/Q); % 生成透镜的透过率函数

lens=pad2center(lens,M,M); % 对透镜进行零填充至 M*M

% 计算从透镜前表面到透镜后表面的透射
u2=u1.*lens;

%%
z2=1/(1/lens_F-1/z1); % 利用成像公式计算像距
z2=z2+2; % 离焦 2 mm

N=round((lambda*z2)/id^2); % 第2次衍射计算时，计算窗口内的采样点数

% 按照计算窗口采样点数，对透镜后表面的光场进行裁剪或零填充
if N<=M
    u2=cutcenter(u2,N,N);
else
    u2=pad2center(u2,N,N);
end

% 计算从透镜后表面到像面的衍射
u3=fresnelas(u2); % 利用角谱法计算
% u3=fresnelsft(u2); % 利用单次傅里叶变换算法计算

%%
S=ceil(R*z2/z1); % 按照放大倍率计算像的行列数
ui=cutcenter(u3,S,S); % 对计算结果进行裁剪

ui_abs=abs(ui);
ui_abs=rot90(ui_abs,2); % 将像旋转180度以方便对比

figure,
subplot(1,2,1),imshow(forg);
subplot(1,2,2),imshow(ui_abs);

