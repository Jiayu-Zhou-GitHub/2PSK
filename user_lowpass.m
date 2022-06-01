function [b,a] = user_lowpass(fs,f_symbol)
%USER_LOWPASS 此处显示有关此函数的摘要
%   此处显示详细说明
Ap=2;As=15;%设置通带允许最大衰减设置为2dB，阻带应达到的最小衰减为15dB
Nn=1000;%抽样次数
F=f_symbol+0.5*f_symbol;%当前频率
F_sh=f_symbol*2;%阻带上限截止频率
%MATLAB工具函数常采用标准化频率，wp和ws的取值范围为0～1
wp=F*2/fs; ws=F_sh*2/fs;%所以通带截止频率为wp/pi ，阻带截止频率为ws/pi 
[N,wn]=buttord(wp,ws,Ap,As);%计算巴特沃斯滤波器阶次和截止频率
[b,a]=butter(N,wn);%频率变换法设计巴特沃斯低通滤波器
end

