clear all
%% BULLISH CALL SPREAD PLOT
fileID = fopen('../data/pricesCall.txt', 'r');
fileID1 = fopen('../data/bsbAskCall.txt', 'r');
fileID2 = fopen('../data/bsbBidCall.txt', 'r');
fileID3 = fopen('../data/bsAskCall.txt', 'r');
fileID4 = fopen('../data/bsBidCall.txt', 'r');
fileID5 = fopen('../data/bsMidCall.txt', 'r');

formatSpec = '%f';
prices = fscanf(fileID, formatSpec);
bsbAsk = fscanf(fileID1, formatSpec);
bsbBid = fscanf(fileID2, formatSpec);
bsAsk = fscanf(fileID3, formatSpec);
bsBid = fscanf(fileID4, formatSpec);
bsMid = fscanf(fileID5, formatSpec);
figure
plot(prices, bsbAsk, '-k')
hold on
plot(prices, bsbBid, '-k')
hold on
plot(prices, bsAsk, '-.r')
hold on
plot(prices, bsBid, '-.r')
hold on
plot(prices, bsMid, '-.b')
hold on
title('Pricing of a Bull Call Spread')
legend('BSB Ask Price', 'BSB Bid Price', 'BS Ask Price', 'BS Bid Price', 'BS mid vol Price');
xlabel('Underlying Asset Price');
ylabel('Derivative Price');

%% BSB,BS COMPARISON PLOT ON EXTREME VOL

fileID6 = fopen('../data/bsLowCall.txt', 'r');
fileID7 = fopen('../data/bsHighCall.txt', 'r');
bsLow = fscanf(fileID6, formatSpec);
bsHigh = fscanf(fileID7, formatSpec);

figure
plot(prices, bsbAsk, '-k')
hold on
plot(prices, bsbBid, '-k')
hold on
plot(prices, bsLow, '-.r')
hold on
plot(prices, bsHigh, '-.r')
hold on

title('BSB and BS Comparison on extreme volatilities')
legend('BSB Ask Price', 'BSB Bid Price', 'BS min vol Price', 'BS max vol Price');
xlabel('Underlying Asset Price');
ylabel('Derivative Price');


%% CALENDAR SPREAD PLOT
fileID = fopen('../data/pricesClnd.txt', 'r');
fileID1 = fopen('../data/bsbAskClnd.txt', 'r');
fileID2 = fopen('../data/bsbBidClnd.txt', 'r');
fileID3 = fopen('../data/bsAskClnd.txt', 'r');
fileID4 = fopen('../data/bsBidClnd.txt', 'r');
fileID5 = fopen('../data/bsMidClnd.txt', 'r');


prices = fscanf(fileID, formatSpec);
bsbAsk = fscanf(fileID1, formatSpec);
bsbBid = fscanf(fileID2, formatSpec);
bsAsk = fscanf(fileID3, formatSpec);
bsBid = fscanf(fileID4, formatSpec);
bsMid = fscanf(fileID5, formatSpec);
figure
plot(prices, bsbAsk, '-k')
hold on
plot(prices, bsbBid, '-k')
hold on
plot(prices, bsAsk, '-.r')
hold on
plot(prices, bsBid, '-.r')
hold on
plot(prices, bsMid, '-.b')
hold on
title('Pricing of a Calendar Spread')
legend('BSB Ask Price', 'BSB Bid Price', 'BS Ask Price', 'BS Bid Price', 'BS mid vol Price');
xlabel('Underlying Asset Price');
ylabel('Derivative Price');



