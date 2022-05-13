%wjq-2019-12-13-608data
%代码目的：初步测试处理数据，导成mat格式
clc
clear
close all
subfunction_path1='.\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.wav';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
% % //======保存图像至指定文件夹===============  
save_directory = ['608matData',date];  %频谱图存储文件夹

if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
   
tic
for i_file=1:length(fname)
%     Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
%     Data=V2Pa_Universal(Data,kulite_transform_ab);
   [Data,Fs]= audioread(fullfile(location,'Database',char(fname(i_file)))); %选择文件导入数据
   %% 对Data进行重组排序
end
fname{1}
