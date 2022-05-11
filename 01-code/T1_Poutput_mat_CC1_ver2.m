
%wjq-2019-11-18
%代码目的：前处理，准备数据，并存储自相关，互相关，DATA结果
clc
clear
% subfunction_path1='H:\动画预测一体化程序通用版_旋转机匣_ver1\subfunction\subfunction_1';
% addpath(subfunction_path1);
%location='E:\Jiaqi-SJTU-DOIT\Database\实验22-2020-模拟声源-16mic - 副本';
location='实验22-2020-模拟声源-1mic';
%[fname,location]=uigetfile({'*.txt';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'/','参数说明','/','parameter.mat']); %选择文件导入数据
disp(Note);
    L_signal =fs*testPeriod;
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft =2^(ceil(log2(L_seg))+1); 
    
fr=[1000,2000,3400];fname=[];
m=[1:16];
for i_m=1:16
    for i_f=1:3
        for i_fi=1
        %fname{i_m,i_f}{1,i_fi}=strcat('m',int2str(m(i_m)),'_',num2str(fr(i_f)),'Hz_Rotate-No-',int2str(i_fi),'.txt');
        fname{i_m,i_f}{1,i_fi}=strcat('mic',int2str(m(i_m)),'_',num2str(fr(i_f)),'_Rotate-No-',int2str(i_fi),'.txt');        
        end
    end
end
%mic1_1000_Rotate-No-1

for i_m=1:16
    for i_f=1:3
      for i_file=1:length(fname{i_m,i_f})
          close all
          Data = importdata(fullfile(location,char(fname{i_m,i_f}(i_file)))); %选择文件导入数据
          Data=V2Pa_Universal(Data,kulite_transform_ab);
          Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
          DATA{i_file}=Data(:,1:13);
    
        for k=1:12
            [temp,freq] = cpsd(DATA{i_file}(:,k),DATA{i_file}(:,13),Wind,Noverlap,Nfft,fs);
            %[temp1,freq1] = cpsd(DATA{i_file}(:,1:12),DATA{i_file}(:,1:12),Wind,Noverlap,Nfft,fs);
            CC1(:,i_file+length(fname{i_m,i_f})*(k-1)) = temp;
            %CC2{i_file} =temp1;
        end
end

%[Pulse,rotor_speed]=keyRotation(DATA{1}(:,end),fs);
%h=figure     
%plot(freq1,10*log10(sum(cell2mat(CC2),2)/length(fname)/12/4e-10))
%title(['Spetrum-Average:', num2str(rotor_speed)])
%saveas(h,['SpetrumAverage', num2str(rotor_speed),'.fig'])

chr=['CC1_if1_m',num2str(m(i_m)),'_f',num2str(i_f),'.mat']
save(chr,'CC1')
end

end

%输出p1-16；f1-3

function DATA=V2Pa_Universal(Data,kulite_transform_ab)
    %first,check the size of Data and kulite_transform_ab,is same or not
    if size(Data,2)~=size(kulite_transform_ab,1)
        disp('转换参数不正确！！')
    end
     for k=1:size(Data,2)  
         DATA(:,k)= Data(:,k)*kulite_transform_ab(k,1)+kulite_transform_ab(k,2);%传感器B1
     end
       
end

 