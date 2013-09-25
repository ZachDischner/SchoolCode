clc;clear all;close all

mystery=load('mystery.txt');
mysteryn = PRNNorm(mystery);

for sat = 1:37
    PRN(sat,:) = PRNNorm(PRNGenerator(sat));
    [lags(sat,:), R(sat,:)] = cyc_corr2(PRN(sat,:),mysteryn');
    Rmax(sat) = max(R(sat,:));
    if (Rmax(sat) == 1)
        fprintf('Satellite is %d, lag is %d',sat,lags(sat))
    end

end