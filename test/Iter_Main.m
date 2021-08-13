clc;
clear all;

%% Parameters Setting
N=512;
M=512;
mes=0.8;
ADC_switch=1;
TestNum=1e2;     %实验次数
IterNum=20;      %算法迭代次数
modType='QAM';   
mod_size=2;
bit=6;
snr=12;

%% Load Parameters
Input.N=N;
Input.M=M;
Input.ADC_switch=ADC_switch;
Input.IterNum=IterNum;
Input.modType=modType;
Input.mod_size=mod_size;
Input.bit=bit;
Input.mes=mes;
Input.nuw=10^(-snr/10);

GEC_MSE_Error=zeros(TestNum,IterNum);
GEC_MSE_Mean=zeros(1,IterNum);
GECSE_MSE_Error=zeros(TestNum,IterNum);
GECSE_MSE_Mean=zeros(1,IterNum);

for kk=1:TestNum
    obj=MIMO_system(Input);
    GEC_MSE_Error(kk,:)=GECSR(obj,Input);
    %GECSE_MSE_Error(kk,:)=GECSR_SE(obj,Input);
    if mod(kk,TestNum/10)==0
        disp(kk/TestNum*10);
    end
end

for ii=1:IterNum
    GEC_MSE_Mean(ii)=mean(GEC_MSE_Error(:,ii));
    %GECSE_MSE_Mean(ii)=mean(GECSE_MSE_Error(:,ii));
end



iter=1:IterNum;
semilogy(iter,   GEC_MSE_Mean, 'LineStyle', '-','LineWidth', 1,  'Color','b', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   hold on;
%semilogy(iter,   GECSE_MSE_Mean, 'LineStyle', 'none','LineWidth', 1,  'Color','b', 'Marker', '*', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   hold on;
legend('GEC-SR'); hold on;
%legend('GEC-SR','SE'); hold on;

xlabel('Iter');
ylabel('MSE');
