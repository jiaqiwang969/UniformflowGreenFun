
clc
clear  
%close all
mode=[1:16];
for zH=0.3%[0.02 0.1 0.3 0.5 1.5 5]
    for i_f=1:3%:3
       for i_m=1:16
close all    
chr.location='.';
frequency=[1000,2000,3400];df=1/(2.56/2);iFreq=round(frequency/df);
chr.g0=['Gnm0_f',num2str(i_f),'_zH',num2str(zH),'.mat'];Gnm0=importdata(fullfile(chr.location,chr.g0));
chr.g1=['Gnm1_f',num2str(i_f),'_zH',num2str(zH),'.mat'];Gnm1=importdata(fullfile(chr.location,chr.g1));
chr.g2=['Tnm1_f',num2str(i_f),'_zH',num2str(zH),'.mat'];Tnm1=importdata(fullfile(chr.location,chr.g2));
chr.g16=['CC1_if1_m',num2str(mode(i_m)),'_f',num2str(i_f),'.mat'];  CC1=importdata(fullfile(chr.location,chr.g16));
% CC1=zeros(size(importdata(fullfile(chr.location,['CC1_p',num2str(1),'_','f',num2str(i_f),'.mat']))));
% for ip=1:16
% CC1=CC1+importdata(fullfile(chr.location,['CC1_p',num2str(ip),'_','f',num2str(i_f),'.mat']));
% end
%CC1=importdata(fullfile(chr.location,['CC1_p',num2str(11),'_','f',num2str(3),'.mat']))+importdata(fullfile(chr.location,['CC1_p',num2str(1),'_','f',num2str(2),'.mat']));


%CC1=CC1_1+CC1_16;
%% Parameter
zR = 0;                                 %重建距离
Ratio=0.6; 
rT=0.185;
%% ----------重建面------------------------------------------
i_r=30;i_cirshift=25;i_equal=16;i_theta=i_equal*i_cirshift;%16个等间隔传感器，因此为16的整数倍，方便累加，无需差值
theta2=(0:i_theta-1)*(360/i_theta);
xr2=linspace(Ratio*rT,rT,i_r);%不沾边，防止Gnm为NA
xr=xr2'*cos(theta2*pi/180);yr=xr2'*sin(theta2*pi/180);%Xnet
zr=0*repmat(1,size(xr));
%%
%CSM is by 16 equal sound


CSM_1=Gnm1(614,:).'*conj(Gnm1(614,:)); 
CSM_2=Tnm1(614,:).'*conj(Tnm1(614,:)); 


%%
CSM_exp=CC1(iFreq(i_f),:).'*conj(CC1(iFreq(i_f),:));%互功率谱矩阵Mo-(2.37）
%Spp_G=CSM_1;Spp_T=CSM_2;Spp_exp=CSM_exp;

%Spp_G=(CSM_1 + CSM_1')./2;Spp_T=(CSM_2 + CSM_2')./2;
%Spp_exp=(CSM_exp + CSM_exp')./2;

Gnm{1}=Gnm0;
Gnm{2}=Gnm1;
Gnm{3}=Tnm1;
toolName={'Gnm0';'Gnm1';'Tnm1'};
sourceName={'Spp_G','Spp_T','Spp_e_x_p'};
viewDir=[0 -90;0 -90;0 -90;180 90;180 90;180 90;180 90;180 90;180 90;];


for i_k=1:3
% H_discretebeam=Gnm{i_k}';
% numberOfeqSources=i_r*i_theta;
% q_sources = size(numberOfeqSources,1);
% steering_matrix = zeros(size(H_discretebeam));
% for index_sources = 1:numberOfeqSources
%      % direct steering vector is vector length of M
%      norm_vector = norm(H_discretebeam(:,index_sources),2);
%      % steering_vector = H_discretebeam(:,index_sources); % H_discretebeam: 30 * 483 
%      steering_vector = H_discretebeam(:,index_sources)/(norm_vector);
%      steering_matrix(:,index_sources) = steering_vector;
%      q_beamforming_G(index_sources) =  steering_vector'*Spp_G*steering_vector;
%      q_beamforming_T(index_sources) =  steering_vector'*Spp_T*steering_vector;     
%      q_beamforming_exp(index_sources) =  steering_vector'*Spp_exp*steering_vector;
% 
% end

%CSM(:,:,k)=Gnm{k}(80*i_r+kk,:).'*conj(Gnm{k}(80*i_r+kk,:));
%CSM_diag0(:,:,k)=CSM(:,:,k);%-diag(diag(CSM(:,:,k)));%剔除噪声贡献明显的自功率谱成分
% Yn=sum(bsxfun(@times,conj(Gnm{knm}),permute(CSM*Gnm{knm}.',[2,1])),2);%消除-1的影响
% q_recon_beamforming2 = reshape(abs(Yn),i_r,i_theta);% sum for each row

q_beamforming_G=sum(bsxfun(@times,conj(Gnm{i_k}),permute(CSM_1*Gnm{i_k}.',[2,1])),2);%消除-1的影响
q_beamforming_T=sum(bsxfun(@times,conj(Gnm{i_k}),permute(CSM_2*Gnm{i_k}.',[2,1])),2);%消除-1的影响
q_beamforming_exp=sum(bsxfun(@times,conj(Gnm{i_k}),permute(CSM_exp*Gnm{i_k}.',[2,1])),2);%消除-1的影响


% q_beamforming_G=sum(bsxfun(@times,conj(Gnm{i_k}.^(-1)),permute(CSM_1'*(Gnm{i_k}.^(-1)).',[2,1])),2).';
% q_beamforming_T=sum(bsxfun(@times,conj(Gnm{i_k}.^(-1)),permute(CSM_2'*(Gnm{i_k}.^(-1)).',[2,1])),2).';
% q_beamforming_exp=sum(bsxfun(@times,conj(Gnm{i_k}.^(-1)),permute(CSM_exp'*(Gnm{i_k}.^(-1)).',[2,1])),2).';
xr=xr2'*cos(theta2*pi/180);yr=xr2'*sin(theta2*pi/180);%Xnet

location=find(ismissing(q_beamforming_G)==1);
q_beamforming_G(location)=1e-10;
q_recon_beamforming_G = reshape(abs(q_beamforming_G),i_r,i_theta);% sum for each row
h1(i_k)=figure;%view(viewDir(i_k,:));
%Qref = max(max(q_recon_beamforming));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
%h1 = pcolor(xr,yr,10*log10(q_recon_beamforming./Qref));% plot in dB unit,hi is the 
surface(yr,xr,q_recon_beamforming_G,'LineStyle','none');% plot in dB unit,hi is the 
title(['Simu:Spp_G-by ',toolName{i_k}])

location=find(ismissing(q_beamforming_T)==1);
q_beamforming_T(location)=1e-10;
q_recon_beamforming_T = reshape(abs(q_beamforming_T),i_r,i_theta);% sum for each row
h2(i_k)=figure;%view(viewDir(i_k,:));
%Qref = max(max(q_recon_beamforming));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
%h1 = pcolor(xr,yr,10*log10(q_recon_beamforming./Qref));% plot in dB unit,hi is the 
surface(yr,xr,q_recon_beamforming_T,'LineStyle','none');% plot in dB unit,hi is the 
title(['Simu:Spp_T-by ',toolName{i_k}])

location=find(ismissing(q_beamforming_exp)==1);
q_beamforming_exp(location)=1e-10;
q_recon_beamforming_exp = reshape(abs(q_beamforming_exp),i_r,i_theta);% sum for each row
h3(i_k)=figure;%view(viewDir(i_k,:));
%Qref = max(max(q_recon_beamforming));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
%h1 = pcolor(xr,yr,10*log10(q_recon_beamforming./Qref));% plot in dB unit,hi is the 
surface(yr,xr,q_recon_beamforming_exp,'LineStyle','none');% plot in dB unit,hi is the 
title(['Exp:Spp_exp-by ',toolName{i_k}])


end
ax(1) = get(h1(1), 'CurrentAxes'); ax(4) = get(h1(2), 'CurrentAxes'); ax(7) = get(h1(3), 'CurrentAxes');                                               
ax(2) = get(h2(1), 'CurrentAxes'); ax(5) = get(h2(2), 'CurrentAxes'); ax(8) = get(h2(3), 'CurrentAxes');   
ax(3) = get(h3(1), 'CurrentAxes'); ax(6) = get(h3(2), 'CurrentAxes'); ax(9) = get(h3(3), 'CurrentAxes');

H=figure;set(gcf,'outerposition',get(0,'screensize'));
suptitle(['f',num2str(i_f),'p16im',num2str(i_m),'zH',num2str(zH)])
for iloop = 1:9
    subp=subplot(3,3,iloop)                                                          % 子图循环
    axChildren = get(ax(iloop),'Children');                                     % 获取axes所有子对象
    copyobj(axChildren, gca); view(subp,viewDir(iloop,:));           % 复制对象到子图的axes 
    title(['Simu:',sourceName{mod(iloop-1,3)+1},'-by ',toolName{ceil(iloop/3)}])
end
saveas(H,['fig\','f',num2str(i_f),'if1im',num2str(i_m),'zH',num2str(zH),'.fig'])
saveas(H,['fig\','f',num2str(i_f),'if1im',num2str(i_m),'zH',num2str(zH),'.png'])
    end
    end
end