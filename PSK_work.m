clc;
clear;
close all;

%%示波器系统
fs=2000000;                         %采样频率
display_t=0.004;                    %数据观测时间

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%调制过程%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name","调制过程");
%%生成单极性NRZ基带信号
f_symbol=25000;                    %码元频率
N=f_symbol*display_t;           %在数据观测时间内的码元个数
T=1/f_symbol;                   %码元持续时间

N_sample=fs*T;                  %每个码元内的采样点数
dt=T/N_sample;                  %时域的采样间隔
a=randi(2,1,N)-1;               %产生单极性的数字随机序列
NRZ=zeros(1,N*N_sample);        %声明变量的空间
%通过for循环将数字序列变为时域的波形
for i=1:N
    for k=1:N_sample
        NRZ((i-1)*N_sample+k)=a(i);
    end
end
t=0:dt:N*N_sample*dt-dt;
subplot(4,1,1);
plot(t,NRZ);
title("单极性NRZ基带信号");
axis([0 display_t -0.5 1.5]);
xlabel('s/t');
ylabel('幅值');

%%生成载波信号
fc=125000;       %载波频率
wc=2*pi*fc;  %载波角频率

dt=1/fs;    %采样间隔
T=1;        %载波观测时间
wt=0:dt:display_t-dt;%与NRZ信号等长
carrier=sin(wc*wt);
subplot(4,1,2);
plot(wt,carrier);
title("载波信号");
axis([0 display_t -1.5 1.5]);
xlabel('s/t');
ylabel('幅值');

%%2PSK调制
PSK_s=zeros(1,display_t*fs);        %声明变量的空间
for t=1:display_t*fs
    if NRZ(t)==1
        PSK_s(t)=sin(wc*t/fs);
    elseif NRZ(t)==0
        PSK_s(t)=sin(wc*t/fs+pi);  
    end
end
subplot(4,1,3);
plot(wt,PSK_s);
title("2PSK信号");
axis([0 display_t -1.5 1.5]);
xlabel('s/t');
ylabel('幅值');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%加入高斯白噪声%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snr=15;%信噪比
px_dBW=0;%噪声功率
PSK=awgn(PSK_s,snr,px_dBW);%加入噪声
subplot(4,1,4);
plot(wt,PSK);
title("加入噪声的2PSK信号");
axis([0 display_t min(PSK)-0.5 max(PSK)+0.5]);
xlabel('s/t');
ylabel('幅值');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%功率谱%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name","功率谱");
%%基带信号功率谱
window=boxcar(length(NRZ)); %矩形窗
[NRZ_Pxx,f]=periodogram(NRZ,window,2048,fs); %直接法
subplot(2,1,1);
plot(f,10*log10(NRZ_Pxx));
title("基带信号功率谱");
axis([0 fs/2 min(10*log10(NRZ_Pxx)) max(10*log10(NRZ_Pxx))]);

%%已调信号功率谱
window=boxcar(length(PSK)); %矩形窗
[PSK_Pxx,f]=periodogram(PSK,window,2048,fs); %直接法
subplot(2,1,2);
plot(f,10*log10(PSK_Pxx));
title("已调信号功率谱");
axis([0 fs/2 min(10*log10(PSK_Pxx)) max(10*log10(PSK_Pxx))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%带通滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%巴特沃斯带通滤波器
figure("Name","带通滤波器设计");
[b,a]=user_bandpass(fs,fc,f_symbol);%带通滤波器
[db,mag,pha,grd,w]=freqz_m(b,a); %求出巴特沃斯带通滤波器幅频、幅度、相位等
subplot(2,1,1);
plot(w*fs/(2*pi),mag);
xlabel('f/Hz');%频率（HZ）
ylabel('幅度/dB');
axis([0,2*fc,0,1.5]);%设置标尺范围
title('数字带通巴特沃斯滤波器')
subplot(2,1,2);
plot(w*fs/(2*pi),180/pi*unwrap(pha));
xlabel('f/Hz');
ylabel('相位/角度');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%带通滤波观察%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name","带通滤波观察");
%2PSK时域
subplot(2,2,1);
plot(wt,PSK);
title("2PSK时域");
axis([0 display_t min(PSK)-0.5 max(PSK)+0.5]);
xlabel('s/t');
ylabel('幅值');
%2PSK频域
f_PSK=fft(PSK,length(PSK));
ff=linspace(0,fs/2,length(PSK)/2);
subplot(2,2,2);
plot(ff,abs(f_PSK(1:length(PSK)/2))/(length(PSK)/2));
title("2PSK频域");
%滤波
o_PSK=filter(b,a,PSK);
%2PSK经过滤波器后时域
subplot(2,2,3);
plot(wt,o_PSK);
title("2PSK经过带通滤波器时域");
axis([0 display_t min(o_PSK)-0.5 max(o_PSK)+0.5]);
xlabel('s/t');
ylabel('幅值');
%2PSK经过滤波器后频域
f_PSK=fft(o_PSK,length(o_PSK));
ff=linspace(0,fs/2,length(o_PSK)/2);
subplot(2,2,4);
plot(ff,abs(f_PSK(1:length(o_PSK)/2))/(length(o_PSK)/2));
title("2PSK经过带通滤波器频域");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%与载波相乘%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_PSK=o_PSK.*2.*carrier;

%%%%%%%%%%%%%%%%%%%%%%%%%%%低通滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name","低通滤波设计")
[b,a] = user_lowpass(fs,f_symbol);%低通滤波器
[db,mag,pha,grd,w]=freqz_m(b,a);
subplot(2,1,1)
plot(w*fs/(2*pi),mag);%db为分贝，mag为增益
xlabel('f/Hz');%频率（HZ）
ylabel('幅度');
title('数字低通巴特沃斯滤波器')
subplot(2,1,2);
plot(w*fs/(2*pi),180/pi*unwrap(pha));
xlabel('f/Hz');
ylabel('相位/角度');

%%%%%%%%%%%%%%%%%%%%%%%%%%%低通滤波器观察%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name","低通滤波观察");
%2PSK经过滤波器后时域
subplot(2,2,1);
plot(wt,o_PSK);
title("2PSK经过带通滤波器时域后与载波相乘时域");
axis([0 display_t min(o_PSK)-0.5 max(o_PSK)+0.5]);
xlabel('s/t');
ylabel('幅值');
%2PSK经过滤波器后频域
f_PSK=fft(o_PSK,length(o_PSK));
ff=linspace(0,fs/2,length(o_PSK)/2);
subplot(2,2,2);
plot(ff,abs(f_PSK(1:length(o_PSK)/2))/(length(o_PSK)/2));
title("2PSK经过带通滤波器时域后与载波相乘频域");
%低通滤波
o_NRZ=filter(b,a,o_PSK)*2;

%2PSK经过滤波器后时域
subplot(2,2,3);
plot(wt,o_NRZ);
title("最终通过低通滤波器时域");
axis([0 display_t min(o_NRZ)-0.5 max(o_NRZ)+0.5]);
xlabel('s/t');
ylabel('幅值');
%2PSK经过滤波器后频域
f_PSK=fft(o_NRZ,length(o_NRZ));
ff=linspace(0,fs/2,length(o_NRZ)/2);
subplot(2,2,4);
plot(ff,abs(f_PSK(1:length(o_NRZ)/2))/(length(o_NRZ)/2));
title("最终通过低通滤波器频域");

figure
%2PSK经过滤波器后时域
plot(wt,o_NRZ);
title("最终通过低通滤波器时域");
axis([0 display_t min(o_NRZ)-0.5 max(o_NRZ)+0.5]);
xlabel('s/t');
ylabel('幅值');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%抽样判决%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure("Name","已解调信号与NRZ对比");
subplot(2,1,1);
plot(wt,NRZ);
title("单极性NRZ基带信号");
axis([0 display_t -0.5 1.5]);
xlabel('s/t');
ylabel('幅值');
output = decision(o_NRZ);%抽样判决
subplot(2,1,2);
plot(wt,output);
title("已解调信号");
axis([0 display_t -0.5 1.5]);
xlabel('s/t');
ylabel('幅值');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%眼图绘制%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=1/f_symbol*fs;%计算出一个码元周期内的采样点数
eyediagram(PSK,n);
title("2PSK信号");
eyediagram(o_NRZ,n*2);
title("解调信号");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%误码率曲线绘制%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt=82;
erro_NRZ=zeros(1,length(NRZ)-nt);
num_1=0;
num_0=0;
for n=1:length(NRZ)-nt          %循环做移位处理
    erro_NRZ(n)=NRZ(n);
    if NRZ(n) == 1
       num_1=num_1+1;
    elseif NRZ(n) == 0
        num_0=num_0+1;
    end
end
p1_NRZ=num_1/(num_1+num_0);%先验概率
p0_NRZ=num_0/(num_1+num_0);

erro_out_NRZ=zeros(1,length(NRZ)-nt);
SNR_min=-6;
SNR_max=18;
d_SNR=0.05;
pe_simulation=zeros(1,(SNR_max-SNR_min)/d_SNR+1);
theoretical=zeros(1,(SNR_max-SNR_min)/d_SNR+1);
for SNR=SNR_min:d_SNR:SNR_max
    eorr_PSK=awgn(PSK_s,SNR,px_dBW);%加入噪声
    [b,a]=user_bandpass(fs,fc,f_symbol);%带通滤波器
    erro_o_PSK=filter(b,a,eorr_PSK);%通过带通滤波
    erro_o_PSK=erro_o_PSK.*2.*carrier;%与载波相乘
    [b,a] = user_lowpass(fs,f_symbol);%低通滤波器
    erro_o_NRZ=filter(b,a,erro_o_PSK)*2;%通过低通滤波器
    erro_output = decision(erro_o_NRZ);%抽样判决
    
    num_01=0;
    num_10=0;
    for n=1:length(NRZ)-nt         %把延时去掉（移位处理）
        erro_out_NRZ(n)=erro_output(n+nt);
        if (erro_out_NRZ(n) == 0)&&(erro_NRZ(n)==1)%发1收到0
            num_01=num_01+1;
        elseif (erro_out_NRZ(n) == 1)&&(erro_NRZ(n)==0)%发0收到1
            num_10=num_10+1;
        end 
    end
    
    p01_erro=num_01/num_1;%后验概率,发1收到0
    p10_erro=num_10/num_0;%后验概率,发0收到1
    pe_simulation(uint16((SNR-SNR_min)/d_SNR+1))=p1_NRZ*p01_erro+p0_NRZ*p10_erro;
    theoretical(uint16((SNR-SNR_min)/d_SNR+1))=1/2*erfc(sqrt(10^(SNR/10)));
end
figure("Name","误码率曲线")
% pe_simulation=10*log(pe_simulation);
% theoretical=10*log(theoretical);
dB=SNR_min:d_SNR:SNR_max;
plot(dB,pe_simulation);
hold;
plot(dB,theoretical);
title("误码率曲线");
legend("仿真误码率曲线","理论误码率曲线");
axis([SNR_min SNR_max min(min(pe_simulation),min(theoretical)) max(max(pe_simulation),max(theoretical))]);
xlabel('r/dB');
ylabel('Pe');