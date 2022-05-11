
%wjq-2019-11-18
%����Ŀ�ģ�ǰ����׼�����ݣ����洢����أ�����أ�DATA���
clc
clear
% subfunction_path1='H:\����Ԥ��һ�廯����ͨ�ð�_��ת��ϻ_ver1\subfunction\subfunction_1';
% addpath(subfunction_path1);
%location='E:\Jiaqi-SJTU-DOIT\Database\ʵ��22-2020-ģ����Դ-16mic - ����';
location='ʵ��22-2020-ģ����Դ-1mic';
%[fname,location]=uigetfile({'*.txt';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'/','����˵��','/','parameter.mat']); %ѡ���ļ���������
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
          Data = importdata(fullfile(location,char(fname{i_m,i_f}(i_file)))); %ѡ���ļ���������
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

%���p1-16��f1-3

function DATA=V2Pa_Universal(Data,kulite_transform_ab)
    %first,check the size of Data and kulite_transform_ab,is same or not
    if size(Data,2)~=size(kulite_transform_ab,1)
        disp('ת����������ȷ����')
    end
     for k=1:size(Data,2)  
         DATA(:,k)= Data(:,k)*kulite_transform_ab(k,1)+kulite_transform_ab(k,2);%������B1
     end
       
end

 