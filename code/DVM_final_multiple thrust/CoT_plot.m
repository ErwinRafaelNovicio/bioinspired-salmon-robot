clc
close all
clear all


%load('data_CoT.mat')
load('CoT_temp.mat')
%cot
frequencies=[1,2,5,7,10,12];
heights=[15,20,22,28,30,32];

figure
surf(heights,frequencies,-CoT')
xlabel('max displacement 3rd active link  [mm]')
ylabel('Freq [Hz]')
title('CoT')
%cot
figure
waterfall(heights,frequencies,-CoT')
xlabel('max displacement 3rd active link  [mm]')
ylabel('Freq [Hz]')
title('CoT')