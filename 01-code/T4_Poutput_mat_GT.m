%��ɢ�������ɳ��͹ܵ����ֺ���

clc
clear
%close all
%bug:Ritio=0.01������

%% �������
save_directory='Beamforming������-duct';
mkdir(save_directory)
subfunction_path1='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\chebfun-master';
addpath(genpath(subfunction_path1));
subfunction_path2='./src';
addpath(subfunction_path2);


%% Case1:����Ƶ�ʷ�Χ�����Case2:Ƶ�ʷ�Χ�������

frequency=[1000,2000,3400];

%%
mode=[-50:50];n_len=7;
N =81;Ratio=0.6; [D,r] = cheb(N,Ratio,1); 
Mx=0.0*ones(N+1,1);%Mx=0.9-r.^2*0.9;%Figure11 (b)

%% Parameter
c = 344;                                 %����
zR = 0;                                  %�ؽ�����
zH = [0.3]                               %��������
zH_model=0.02*[1:40]                     %���潨ģ������������̣�����zH��Ϊ������Ȼ��*exp(sqrt(-1)*k_nm(kk)*��zH_model-zH)����ע��ֻ�����Դͬһ��
%zH_model=[0.2 0.3]                      %���潨ģ������������̣�����zH��Ϊ������Ȼ��*exp(sqrt(-1)*k_nm(kk)*��zH_model-zH)����ע��ֻ�����Դͬһ��

rT=0.185;
%% ----------�ؽ���------------------------------------------
i_r=30;i_cirshift=25;i_equal=16;i_theta=i_equal*i_cirshift;%16���ȼ�������������Ϊ16���������������ۼӣ������ֵ
theta2=(0:i_theta-1)*(360/i_theta);
xr2=linspace(Ratio*rT,rT,i_r);%��մ�ߣ���ֹGnmΪNA
xr=xr2'*cos(theta2*pi/180);yr=xr2'*sin(theta2*pi/180);%Xnet
zr=0*repmat(1,size(xr));
Angle=90;%ż������Դ�ķ���90�ȶ�Ӧʵ��

%% ������circle type������������������������������
% aperture size=0.2m,num of microphones is 360
%theta1=(0:length(Sig.fname)*12-1)*(360/length(Sig.fname)/12);
circle_number=6*12;
theta1=(0:circle_number-1)*(360/circle_number); %�������ε������ؽ���Դ������


xh=rT*cos(theta1*pi/180);yh=rT*sin(theta1*pi/180);Nxh=length(xh);%Xmic
zh=zH*ones(1,length(xh));  
xr1=rT*linspace(Ratio,1,i_r).'*cos(linspace(0,2*pi,i_theta));yr1=rT*linspace(Ratio,1,i_r).'*sin(linspace(0,2*pi,i_theta));%Xnet��������ؽ�����һ�µģ���ʵ�����Ѿ���ģ�����Դ���ϵĽ��
%%

%% ---------�������þ�����ʽƵ�򷽷�����---------------
R=[xr(:) yr(:) zr(:)];
H=[xh(:) yh(:) zh(:)];
Dist=sqrt(R.^2*ones(size(H'))+ones(size(R))*(H').^2-2*R*H');
r_pole=linspace(Ratio,1,i_r);%�����λһ��


for k=1:length(frequency)
%ww����ɨ������,cormatrix������ؾ���
%https://blog.csdn.net/qq_36300268/article/details/88739184 
%���˼·��������ȫ����һ��
Gnm0=exp(-1i*2*pi*bsxfun(@times,Dist,frequency(k))/c)./(4*pi*Dist); %steer vector,(2.30)
%����ϣ������Gnm����ܵ��ı߽�����
%���ͣ�Gnm��ʾ

[initialEigValue,mode_enlarge,cutOffLine,len]=wm2initialEigValue(N,D,r,Ratio,Mx,frequency(k)*2*pi*rT/c,mode,n_len);
%bug:��û������������,�����λ���û���֡�
for ki=1:length(initialEigValue)
    if length(initialEigValue{ki})>=2 
        initialEigValue_left{ki}= initialEigValue{ki}(1:floor(length(initialEigValue{ki})/2));
        initialEigValue_right{ki}= initialEigValue{ki}(end-floor(length(initialEigValue{ki})/2)+1:end);
        mode_enlarge_left{ki}=mode_enlarge{ki}(1:floor(length(initialEigValue{ki})/2));
        mode_enlarge_right{ki}=mode_enlarge{ki}(end-floor(length(initialEigValue{ki})/2)+1:end);    
        len_left{ki}=len{ki}(1:floor(length(initialEigValue{ki})/2));
        len_right{ki}=len{ki}(end-floor(length(initialEigValue{ki})/2)+1:end);      
    else
        initialEigValue_left{ki}= [];
        initialEigValue_right{ki}= [];      
        mode_enlarge_left{ki}=[];
        mode_enlarge_right{ki}=[];    
        len_left{ki}=[];
        len_right{ki}=[];  
    end  
end

EigValue_left=cell2mat(initialEigValue_left.');EigValue_right=cell2mat(initialEigValue_right.');
Mode_enlarge_left=cell2mat(mode_enlarge_left.');Mode_enlarge_right=cell2mat(mode_enlarge_right.');
Len_left=cell2mat(len_left);Len_right=cell2mat(len_right);

%��Ѽ��ͼ����������ֵͼ��ֻ�����ǹ̶�Ƶ�����棬��ͬģ̬�����ԣ�textע�͵�Ϊ��ͬ����ģ̬��
figure;plot3(real(EigValue_left),imag(EigValue_left),Mode_enlarge_left,'.r');hold on;plot3(real(EigValue_right),imag(EigValue_right),Mode_enlarge_right,'.b');
view([0 0]);
text(real(EigValue_left),imag(EigValue_left),Mode_enlarge_left,num2str(Len_left'));
text(real(EigValue_right),imag(EigValue_right),Mode_enlarge_left,num2str(Len_right'));
grid on;
para.TPP1=fill3([cutOffLine cutOffLine cutOffLine  cutOffLine],[-20 -20 20 20],[-20 20 20 -20],'r','FaceAlpha',0.1); %�����ɫ        % handle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set(para.TPP1,'Visible','on')  
para.TPP2=fill3([20 -20 -20 20],[0 0 0 0],[-20 -20 20  20],'g','FaceAlpha',0.1); %�����ɫ        % handle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set(para.TPP2,'Visible','on') 
view([16 21]);
%% ������ 
Tm=[];Tv=[];Gm=[];Gv=[];
tic
[Tm,Tv,Gm,Gv]=greenfun_bf_initial(N,D,r,Ratio,Mx,frequency(k)*2*pi*rT/c,Nxh,i_r,i_theta,Mode_enlarge_right,EigValue_right,zH/rT,r_pole,Angle);% ż������Դ
toc
[Gnm1,Tnm1]=greenfun_bf_wall(Tm,Gm,Nxh,i_r,i_theta,Mode_enlarge_left,EigValue_left,r_pole);


chr0=['Gnm0_f',num2str(k),'_zH',num2str(zH),'.mat'];save(chr0,'Gnm0');
chr1=['Gnm1_f',num2str(k),'_zH',num2str(zH),'.mat'];save(chr1,'Gnm1');
chr2=['Tnm1_f',num2str(k),'_zH',num2str(zH),'.mat'];save(chr2,'Tnm1');
chr3=['Gv_f',num2str(k),'_zH',num2str(zH),'.mat'];save(chr3,'Gv');
chr4=['Tv_f',num2str(k),'_zH',num2str(zH),'.mat'];save(chr4,'Tv');


end  
        
        




