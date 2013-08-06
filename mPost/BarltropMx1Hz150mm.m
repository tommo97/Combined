%close all
%clear all
clc
load Barltrop1Hz150mmMy.mat
[a b] = sort(data(:,1));
data_srt = data(b,1:2);
[a b] = unique(data_srt(:,1));
data_srt = data_srt(b,:);
data_srt_int(:,1) = linspace(min(data_srt(:,1)),max(data_srt(:,1)),round((max(data_srt(:,1))-min(data_srt(:,1)))*1000));
data_srt_int(:,2) = interp1(data_srt(:,1),data_srt(:,2),data_srt_int(:,1));
%scatter(data_srt_int(:,1),smooth(data_srt_int(:,2)),'x')