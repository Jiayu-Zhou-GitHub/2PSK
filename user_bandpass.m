function [b,a] = user_bandpass(fs,fc,f_symbol)
%USER_HIGHPASS 此处显示有关此函数的摘要
%   此处显示详细说明
wp_d=fc-f_symbol-8000;%通带下限截止频率
wp_s=fc+f_symbol+8000;%通带上限截止频率

ws_d=fc-f_symbol*2;%阻带下限截止频率
ws_s=fc+f_symbol*2;%阻带上限截止频率
wp=[wp_d wp_s]*2/fs; ws=[ws_d ws_s]*2/fs;%归一化
Ap=3;As=20;%设置通带允许最大衰减设置为3dB，阻带应达到的最小衰减为18dB
[N,wn]=buttord(wp,ws,Ap,As);%计算巴特沃斯滤波器阶次和截止频率
%N为滤波器阶数 wn为滤波器截止频率 数字滤波器
[b,a]=butter(N,wn,'bandpass');%频率变换法设计巴特沃斯带通滤波器
end

