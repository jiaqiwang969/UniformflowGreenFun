%wjq-2019-12-13-608data
%����Ŀ�ģ��������Դ������ݣ�����mat��ʽ
clc
clear
close all
subfunction_path1='.\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.wav';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'\','����˵��','\','parameter.mat']); %ѡ���ļ���������
% % //======����ͼ����ָ���ļ���===============  
save_directory = ['608matData',date];  %Ƶ��ͼ�洢�ļ���

if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('�ļ��д��ڣ�');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
   
tic
for i_file=1:length(fname)
%     Data = importdata(fullfile(location,char(fname(i_file)))); %ѡ���ļ���������
%     Data=V2Pa_Universal(Data,kulite_transform_ab);
   [Data,Fs]= audioread(fullfile(location,'Database',char(fname(i_file)))); %ѡ���ļ���������
   %% ��Data������������
end
fname{1}
