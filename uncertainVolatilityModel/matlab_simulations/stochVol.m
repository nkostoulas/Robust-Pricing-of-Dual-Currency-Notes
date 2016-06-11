%clear all

fileID = fopen('../data/stochVolData.txt','r');
%fileID = fopen('fig4_100000_data.txt', 'r');
formatSpec = '%f';
A = fscanf(fileID, formatSpec);
histogram(A,'Normalization','Probability')  %pdf or probability
title('Probability Distribution of Profit/Loss for BSB Delta Hedging');
%title('Probability Distribution of Profit/Loss for BS Delta Hedging');
xlabel('Profit/Loss');
ylabel('Probability');

%histfit(A)







%produce pdf only from normal fit of histfit
%[mu, sigma] = normfit(A);
%pd = fitdist(A,'normal');
%x_values = -5.8:0.01:5.8;
%y = pdf(pd,x_values);
%figure
%plot(x_values,y,'LineWidth',2)

% 2
% [f,x]=hist(A, 50);%# create histogram from a normal distribution.
% g = normpdf(x, mu, sigma);
% 
% %#METHOD 1: DIVIDE BY SUM
% figure(1)
% bar(x,f/sum(f));hold on
% plot(x,g,'r');hold off
% 
% %#METHOD 2: DIVIDE BY AREA
% figure(2)
% bar(x,f/trapz(x,f));hold on
% plot(x,g,'r');hold off

% 3
% numOfBins = 50;
% [histFreq, histXout] = hist(A, numOfBins);
% figure;
% bar(histXout, histFreq/sum(histFreq));
% xlabel('x');
% ylabel('Frequency (percent)');