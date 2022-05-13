%wjq-2019-12-13-608data
%代码目的：初步测试处理数据，导成mat格式
clc
clear
close all
% subfunction_path1='.\subfunction_1';
% addpath(subfunction_path1);
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

%  Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
%  Data=V2Pa_Universal(Data,kulite_transform_ab);
   %% 对Data进行重组排序
[data,Fs]= audioread(fullfile(location,char(fname(i_file))),'native'); %选择文件导入数据
Data = data(:,mic_transform_ab(:,2));
    DATA=(double(Data)/(2^15)*10)/0.05;       %16位采集卡，测量范围正负10V，分成2^16个格，初始的数据表示范围-2^15~2^15
     figure;plot(DATA)                                                                           %-2^15~2^15整数转为正负10V电压，再除以50mv/pa的灵敏度得到单位为pa的数据
%    waveletPlot_universal(DATA,Fs,2.56)
%    DATA=Data;
   Phase=circle1/180*pi;
%    Phase=[1:32]*360/30/180*pi;;
%   20*log10(max(Data/0.05)/2e-5);
%   max(Data);
%    DATA=DATA./rms(DATA);%对信号做正则化处理,正则化稀疏度更好

   %figure;plot(Data(:,1))
   %[the_freq,freq_dB]=frequencyDomainPlot_dB_no_deal(DATA,Fs,2.56);
   [the_freq1,freq_dB1]=frequencyDomainPlot_dB(DATA,Fs,2.56);   [the_freq1,freq_dB1]=frequencyDomainPlot_dB(DATA,Fs,2.56);

%   max(freq_dB1)
    nk=32;
    
    BPF=the_freq1(find(freq_dB1(:,1)==max(freq_dB1(:,1))));rotor_speed=BPF/17;
    listName={'1BPF';'2BPF';};    mode=[-50:50];
    [GAMMA,freq]=wavemode_calculation_CPSD_phase(DATA,Fs,nk,mode,Phase,rotor_speed*60,save_directory,i_file,fname,1);
    [GAMMA,freq]=wavemode_calculation_CPSD_phase(DATA,Fs,nk,mode,Phase,rotor_speed*60,save_directory,i_file,fname,2);
%     df=freq(2)-freq(1);list=[round(16.5*rotor_speed/df) round(17.5*rotor_speed/df);round(33.5*rotor_speed/df) round(34.5*rotor_speed/df);];
% 
%     for k=1:2
%     [x(k),y(k)]=find(GAMMA==max(max(GAMMA(list(k,1):list(k,2),:))),1);Zdata=GAMMA(x(k),:);Zdata(y(k))=0;[y_1(k)]=find(Zdata==max(Zdata),1);
%     o_RI(k)=x(k)*df;m_RI(k)=y(k)-(mode(end)+1);m_RI_1(k)=y_1(k)-(mode(end)+1);
%     text(m_RI(k),o_RI(k),{[num2str(round(m_RI(k)*10)/10),',',num2str(round(GAMMA(x(k),y(k))))];[num2str(o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60))]});text(m_RI_1(k),o_RI(k),[num2str(m_RI_1(k)),',',num2str(round(GAMMA(x(k),y_1(k))))]);   
%     end
%     for k=1:2
%     h1{k}=figure;bar(mode,GAMMA(x(k),:));hold on
% %     ylim([40 70]);
%     title([char(fname(i_file)),'-',num2str(round(rotor_speed*60)),'rpm','-',listName{k}],'FontSize',14)
%     saveas(h1{k},[save_directory,'\',strrep(char(fname{i_file}),'.wav','-'),'modeplot-',listName{k},'.png'])
%     saveas(h1{k},[save_directory,'\',strrep(char(fname{i_file}),'.wav','-'),'modeplot-',listName{k},'.fig'])

%     end
end




