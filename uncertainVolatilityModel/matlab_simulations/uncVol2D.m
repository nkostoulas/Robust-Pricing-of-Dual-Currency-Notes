clear all
fileID1 = fopen('../data/2Dbsb_spread_L.txt', 'r');
fileID2 = fopen('../data/2Dbsb_spread_U.txt', 'r');
fileID3 = fopen('../data/2Dbs_spread_L.txt', 'r');
fileID4 = fopen('../data/2Dbs_spread_U.txt', 'r');
fileID = fopen('../data/2Dbsb_spread_prices.txt', 'r');

formatSpec = '%f';
prices = fscanf(fileID, formatSpec);
lower_vec = fscanf(fileID1, formatSpec);
upper_vec = fscanf(fileID2, formatSpec);
lower_vec2 = fscanf(fileID3, formatSpec);
upper_vec2 = fscanf(fileID4, formatSpec);

N = size(prices,1);

lower = vec2mat(lower_vec, N);
upper = vec2mat(upper_vec, N);
lower2 = vec2mat(lower_vec2, N);
upper2 = vec2mat(upper_vec2, N);

figure;hold on
colormap([0 0 0;0 0 0;1 0 0;1 0 0])
mesh(lower, ones(N)); % or surf
mesh(upper, ones(N)+1); % or surf
mesh(lower2, ones(N)+2); % or surf
mesh(upper2, ones(N)+3); % or surf
title('Pricing of a Spread of Options on 2 assets');
legend('BSB Bid', 'BSB Ask', 'BS Bid', 'BS Ask');
xlabel('Asset 1 Price');
ylabel('Asset 2 Price');
zlabel('Derivative Price');
view(0,0);
