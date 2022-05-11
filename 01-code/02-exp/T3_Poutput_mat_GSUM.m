%T3 Output GV_SUM
clc
clear
Ratio=0.6; 
zH=0.3;zR = 0;                                 %重建距离
zH_model=0.02*[1:40]                     %仿真建模（代码操作过程：先算zH作为基础，然后*exp(sqrt(-1)*k_nm(kk)*（zH_model-zH)）；注意只针对声源同一侧

for i_f=3
chr.location='./mat';
frequency=[1000,2000,3400];df=1/(2.56/2);iFreq=round(frequency/df);
chr.g0=['Gv_f',num2str(i_f),'_zH',num2str(zH),'.mat'];Gv=importdata(fullfile(chr.location,chr.g0));
chr.g1=['Tv_f',num2str(i_f),'_zH',num2str(zH),'.mat'];Tv=importdata(fullfile(chr.location,chr.g1));
chr.e0=['EigValue_right_zH_-',num2str(zH),'_f',num2str(i_f),'.mat'];EigValue_right=importdata(fullfile(chr.location,chr.e0));
end



rT=0.185;
%% ----------重建面------------------------------------------
i_r=30;i_cirshift=25;i_equal=16;i_theta=i_equal*i_cirshift;%16个等间隔传感器，因此为16的整数倍，方便累加，无需差值
theta2=(0:i_theta-1)*(360/i_theta);
xr2=linspace(Ratio*rT,rT,i_r);%不沾边，防止Gnm为NA
xr=xr2'*cos(theta2*pi/180);yr=xr2'*sin(theta2*pi/180);%Xnet
zr=0*repmat(1,size(xr));
Angle=90;%偶极子声源的方向；90度对应实验

%% ———circle type———————————————
% aperture size=0.2m,num of microphones is 360
%theta1=(0:length(Sig.fname)*12-1)*(360/length(Sig.fname)/12);
circle_number=6*12;
theta1=(0:circle_number-1)*(360/circle_number); %测量环形点阵，在重建声源面的左侧


xh=rT*cos(theta1*pi/180);yh=rT*sin(theta1*pi/180);Nxh=length(xh);%Xmic
zh=zH*ones(1,length(xh));  
xr1=rT*linspace(Ratio,1,i_r).'*cos(linspace(0,2*pi,i_theta));yr1=rT*linspace(Ratio,1,i_r).'*sin(linspace(0,2*pi,i_theta));%Xnet，测点数重建面是一致的，但实际上已经是模拟的声源面上的结果
%%

%% ---------以下是用矩阵形式频域方法运算---------------
R=[xr(:) yr(:) zr(:)];
H=[xh(:) yh(:) zh(:)];
Dist=sqrt(R.^2*ones(size(H'))+ones(size(R))*(H').^2-2*R*H');
r_pole=linspace(Ratio,1,i_r);%脉冲点位一列


%%

tic
[Gv_sum,Tv_sum]=greenfun_bf_ref(Tv,Gv,zH,zH_model/rT,Nxh,i_r,i_theta,EigValue_right,r_pole);
toc
%% 集成16个等间隔分布，偶极子发声
Gv_sum_eq=[];Tv_sum_eq=[];GV_SUM=[];TV_SUM=[];
%Theta=[1:16]
for kz=1:length(zH_model)
    for kr=1:i_r
             GV_SUM{kz}{kr,1}=zeros(size(Gv_sum{kz}{1, 1}));
             TV_SUM{kz}{kr,1}=zeros(size(Tv_sum{kz}{1, 1}));
        for ke=1:i_equal
          Gv_sum_eq{ke,kz}{kr}=circshift(Gv_sum{kz}{kr},(ke-1)*i_cirshift,2)*exp(-sqrt(-1)*7*(ke-1)*2*pi/i_equal);
          Tv_sum_eq{ke,kz}{kr}=circshift(Tv_sum{kz}{kr},(ke-1)*i_cirshift,2)*exp(-sqrt(-1)*7*(ke-1)*2*pi/i_equal);
          GV_SUM{kz}{kr,1}=GV_SUM{kz}{kr,1}+Gv_sum_eq{ke,kz}{kr};
          TV_SUM{kz}{kr,1}=TV_SUM{kz}{kr,1}+Tv_sum_eq{ke,kz}{kr};
        end
    end
end

%%
i_e=1%:i_equal
for N=15
h=MATLAB4geomTurbo(xr,yr,zr,xh,yh,zh,xr1,yr1);

%N=20;
 xishu_a=1600/(max(max(abs(TV_SUM{1,1}{N})))-min(min(abs(TV_SUM{1,1}{N}))));
 xishu_b=1600;
%axes = axes('Parent',h1,'Position',[0.07 0.20 0.40 0.70]);hold(axes1,'on');
 %
for kz=1:length(zH_model)
% figure
surface(yr1*1000,xr1*1000,repmat(-zH_model(kz),size(xr1))*1000-171,imag(TV_SUM{i_e,kz}{N})*xishu_a-xishu_b,'LineStyle','none');%title('Gv');view(axes1,[90 90]);
% shading interp
% alpha(0.8)
%colormap (summer)
end
end