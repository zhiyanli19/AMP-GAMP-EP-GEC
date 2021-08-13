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
bit=8;
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


EP_MSE_Error = zeros(TestNum, IterNum);
EP_MSE_Mean = zeros(1, IterNum);

AMP_MSE_Error = zeros(TestNum, IterNum);
AMP_MSE_Mean = zeros(1, IterNum);

GAMP_MSE_Error = zeros(TestNum, IterNum);
GAMP_MSE_Mean = zeros(1, IterNum);
%BER_MSE_Error = zeros(test_num, 1);
%SER_MSE_Error = zeros(test_num, 1);
j = 1;
for i = 1: TestNum
    obj=MIMO_system(Input);
    %obj = EP_MIMO_newSystem(Input);
    GEC_MSE_Error(i,:)=GECSR(obj,Input);
    [EP_MSE_Error(i, :)] = EP_Iter_new2(obj, Input);
    [EP_MSE_Error(i, :)] = EPnew(obj, Input);
    AMP_MSE_Error(i, :) = AMP_Iter(obj, Input);
    GAMP_MSE_Error(i, :) = GAMP_Iter(obj, Input);
    if mod(i, TestNum / 10) == 0
        disp(i / TestNum * 10);
    end
end
for i = 1: IterNum
    EP_MSE_Mean(i) = mean(EP_MSE_Error(:, i));
    GEC_MSE_Mean(i) = mean(GEC_MSE_Error(:,i));
    AMP_MSE_Mean(i) = mean(AMP_MSE_Error(:,i));
    GAMP_MSE_Mean(i) = mean(GAMP_MSE_Error(:,i));
end

iter=1:IterNum;
semilogy(iter,   GEC_MSE_Mean, 'LineStyle', '-','LineWidth', 1,  'Color','b', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   hold on;
semilogy(iter,   EP_MSE_Mean, 'LineStyle', '-','LineWidth', 1,  'Color','r', 'Marker', '*', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   hold on;
semilogy(iter,   AMP_MSE_Mean, 'LineStyle', '-','LineWidth', 1,  'Color','g', 'Marker', '+', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   hold on;
semilogy(iter,   GAMP_MSE_Mean, 'LineStyle', '-','LineWidth', 1,  'Color','r', 'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   hold on;
%legend('GEC-SR'); hold on;
%legend('GEC','EP'); hold on;
%legend('AMP','GAMP'); hold on;
%legend('GEC-SR','EP','AMP'); hold on;
legend('GEC-SR','EP','AMP','GAMP'); hold on;
xlabel('Iter');
ylabel('MSE');