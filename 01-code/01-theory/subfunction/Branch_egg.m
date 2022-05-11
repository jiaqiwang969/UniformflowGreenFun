
clc
clear


%% 输出配置
save_directory='Beamforming报告结果-duct';
mkdir(save_directory)
subfunction_path1='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\chebfun-master';
addpath(genpath(subfunction_path1));

frequency=[100:100:2900]


%% Parameter
c = 344;                              %声速
zR = 0;                                 %重建距离
zH = [0.02]                                %测量距离
%% ———circle type———————————————
% aperture size=0.2m,num of microphones is 360
%theta1=(0:length(Sig.fname)*12-1)*(360/length(Sig.fname)/12);
theta1=(0:30*12-1)*(360/30/12);

rT=0.2;
xh=rT*cos(theta1*pi/180);yh=rT*sin(theta1*pi/180);Nxh=length(xh);%Xmic
zh=-zH*ones(1,length(xh));

%% ----------重建面------------------------------------------
i_r=30;i_theta=360;%net
theta2=(0:i_theta-1)*(360/i_theta);
xr2=linspace(0.3*rT,rT,i_r);%不沾边，防止Gnm为NA
xr=xr2'*cos(theta2*pi/180);yr=xr2'*sin(theta2*pi/180);%Xnet
zr=0*repmat(1,size(xr));


mode=[-50:50];n_len=7;
N =81;Ratio=0.3; [D,r] = cheb(N,Ratio,1); 
Mx=0.6*ones(N+1,1);
%Mx=0.9-r.^2*0.9;%Figure11 (b)

[initialEigValue,mode_enlarge]=wm2initialEigValue(N,D,r,Ratio,Mx,frequency*2*pi*rT/c,mode,n_len);
for ki=[[1:12] [90:101]] %这里就有可能产生误差，
    initialEigValue{ki}=[];
    mode_enlarge{ki}=[];
end
%bug:还没有区分上下游,上下游还是没区分，这里是认为手动对半分。误差有点大。
for ki=1:length(initialEigValue)
    if length(initialEigValue{ki})>=2 
        initialEigValue_left{ki}= initialEigValue{ki}(1:floor(length(initialEigValue{ki})/2));
        initialEigValue_right{ki}= initialEigValue{ki}(end-floor(length(initialEigValue{ki})/2)+1:end);
        mode_enlarge_left{ki}=mode_enlarge{ki}(1:floor(length(initialEigValue{ki})/2));
        mode_enlarge_right{ki}=mode_enlarge{ki}(end-floor(length(initialEigValue{ki})/2)+1:end);       
    else
        initialEigValue_left{ki}= [];
        initialEigValue_right{ki}= [];      
        mode_enlarge_left{ki}=[];
        mode_enlarge_right{ki}=[];    
    end  
end

EigValue_left=cell2mat(initialEigValue_left.');EigValue_right=cell2mat(initialEigValue_right.');
Mode_enlarge_left=cell2mat(mode_enlarge_left.');Mode_enlarge_right=cell2mat(mode_enlarge_right.');
%画鸭蛋图；还是特征值图，只不过是固定频率下面，不同模态的特性，text注释的为不同周向模态。
figure;plot3(real(EigValue_left),imag(EigValue_left),Mode_enlarge_left,'.');hold on;plot3(real(EigValue_right),imag(EigValue_right),Mode_enlarge_right,'.');
view([0 0]);
text(real(EigValue_left),imag(EigValue_left),Mode_enlarge_left,num2str(Mode_enlarge_left));
text(real(EigValue_right),imag(EigValue_right),Mode_enlarge_left,num2str(Mode_enlarge_right));
      



