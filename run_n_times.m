%% --- clean ---
clear; close all; clc;

%% --- run image_series_segmentation ---
N = inputdlg('Enter the number of experiments:');
N = str2num(N{1});
for i =  1:N
    image_series_segmentation;
    data(:,i) = MediumMI; %save medium stats
    peri{i} = per; %save perimeter
    
%% ---display plots---
figure(2)
plot(MediumMI)
xlabel('Frame');
ylabel('MFI');
hold on
title('medium mean intensity by frame');
end

%add legend
for i=1:N
lgnd{i} = int2str(i);
end
figure(2);
legend(lgnd);

%% --- export to excel ---
cd(selectedPath);
cd .. ;
mkdir('result');
cd('result');
filename = 'data.xlsx';
writematrix(data,'data.xlsx','Sheet',1)

%% --- save segmeneted images ---
for i=1:N
%generate required variables
ImOutFolder=pwd;
ImName= strcat(int2str(i), '.png');
RGBIm=randi(127,100,100,3,'uint8');
LogicalIm=logical(randi([0 1],100,100));
%compose inputs
ImFileOut = fullfile(ImOutFolder, ImName);
ImCell = {RGBIm, LogicalIm};
%create montage and extract color data from the image object
figure(i+2)
h=montage(peri{i});
montage_IM=h.CData;
%write to file
imwrite(montage_IM,ImFileOut);
end