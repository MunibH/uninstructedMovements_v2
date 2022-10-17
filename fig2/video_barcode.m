

% load('C:\Users\munib\Documents\Economo-Lab\data\Video\JEB7\2021-04-29\grayframes_Cam0.mat')
% temp = grayframes_Cam0;

load('C:\Users\munib\Documents\Economo-Lab\data\Video\JEB7\2021-04-29\grayframes_Cam1.mat')
temp = grayframes_Cam1;


%%
close all

figure; 
temp2 = squeeze(mean(temp,1));
% temp2(temp2>120) = 90;

tix = 1:10000;
time = tix ./ 400;
imagesc(time,1:size(temp2,1),temp2(:,tix))