clear all 

%% BULLISH CALL SPREAD PLOT
fileID = fopen('../data/bsbAskDual.txt', 'r');
fileID1 = fopen('../data/bsbBidDual.txt', 'r');
fileID2 = fopen('../data/bsAskDual.txt', 'r');
fileID3 = fopen('../data/bsBidDual.txt', 'r');
fileID4 = fopen('../data/bsMidDual.txt', 'r');
fileID5 = fopen('../data/pricesDual.txt', 'r');
    
formatSpec = '%f';
prices = fscanf(fileID5, formatSpec);
bsbAsk = fscanf(fileID, formatSpec);
bsbBid = fscanf(fileID1, formatSpec);
bsAsk = fscanf(fileID2, formatSpec);
bsBid = fscanf(fileID3, formatSpec);
bsMid = fscanf(fileID4, formatSpec);
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
%title('Pricing of a Bull Call Spread on PRDC notes')
title('Pricing of a Calendar Spread on PRDC notes')
legend('BSB Ask Price', 'BSB Bid Price', 'BS Ask Price', 'BS Bid Price', 'BS mid vol Price');
xlabel('Underlying FX Price');
ylabel('Derivative Price');