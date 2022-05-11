clc
clear
close all

%% Add Subfunction
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));


%%
load('database/mics_loc.mat')
mics_x=mics_loc(:,1).*cos(mics_loc(:,2));
mics_y=mics_loc(:,1).*sin(mics_loc(:,2));
mics_z=mics_loc(:,3);



%% parameters
Angle=0;

f=2000;
a=0.3;
w = f*2*pi/343*0.3;
r1 = 0.2/a;  theta1=0; z1=0/a;   % source term 1
r2 = 0.26/a; theta2=pi; z2=0/a;  % source term 2
%r=1; theta=mics_loc(1:16,2).';z1=mics_loc(:,3).'/a;
r = 0.99; theta = linspace(0,2*pi,16); z1 = [-0.939 -0.839 -0.539]/a;
M = 0;
%% Mode Generator
m = [-50:50];
n = [50];

[Base] = BaseJ1(m,n);

beta=sqrt(1-M^2);
kappa_mn=sqrt(w^2-beta^2*Base.jmn_pm.^2);

Eigm_mn=(-w*M+kappa_mn)/beta^2;  % left running
Eigp_mn=(-w*M-kappa_mn)/beta^2;  % right running


Omagam_mn = w - M*Eigm_mn;   %AIAA-9
Omagap_mn = w - M*Eigp_mn;   %AIAA-9
Qm_mn =  kappa_mn.*(1-m.^2./Base.jmn_pm.^2); %Lowis
Qp_mn =  kappa_mn.*(1-m.^2./Base.jmn_pm.^2);

% Qm_mn =  (Eigm_mn+Omagam_mn*M).*(1-m.^2./Base.jmn_pm.^2); %Rienstra
% Qp_mn =  (Eigm_mn+Omagam_mn*M).*(1-m.^2./Base.jmn_pm.^2);


for km=1:length(m)
    Gmn1(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r1)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qm_mn(:,km))...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigm_mn(:,km)*z1)),n,1,length(z1));
    Gmn2(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r2)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qm_mn(:,km))...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigm_mn(:,km)*z1)),n,1,length(z1));
%     Gmn2(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r1)...
%         ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qp_mn(:,km))...
%         .*besselj(m(km),Base.jmn_pm(:,km)*r)...
%         .*reshape(exp(-i*(Eigp_mn(:,km)*z2)),n,1,length(z2));
end
Gm1=reshape(sum(Gmn1,1),length(r),length(z1),length(m))./(-2*pi*i);       %AIAA-20
Gm2=reshape(sum(Gmn2,1),length(r),length(z1),length(m))./(-2*pi*i);       %AIAA-20
tem1=exp(-1i*m.'*(theta-theta1));
tem2=exp(-1i*m.'*(theta-theta2));

Gw1=reshape(sum(bsxfun(@times,Gm1,reshape(tem1,1,1,length(m),length(theta))),3),length(r),length(z1),length(theta));
Gw2=reshape(sum(bsxfun(@times,Gm2,reshape(tem2,1,1,length(m),length(theta))),3),length(r),length(z1),length(theta));

Gw=reshape(Gw1+Gw2,3,16);



%% Experiments
Data = importdata('database/test20210122_test2_12_2000_1_1.mat');
ref=Data(:,4);
Fs = 16384;            % Sampling frequency
T = 1/Fs;              % Sampling period



data_fft = 2^floor(log2(length(ref)));
data = Data(end+1-data_fft:end,:);
the_freq = [0:data_fft/2.56 - 1]*Fs/data_fft;  %数据频域离散刻度
data_freq = fft(data)*2/data_fft;
data_freq = data_freq(1:data_fft/2.56,:);

% L = length(Data);      % Length of signal
% t = (0:L-1)*T;         % Time vector
% L_signal = length(ref);
% L_seg = round(L_signal/10);
% Wind = hamming(L_seg);
% Noverlap = round(L_seg/2);
% Nfft = 2^(ceil(log2(L_seg))+1);
% for k=1:48
%     [temp,freq] = cpsd(Data(:,k),ref,Wind,Noverlap,Nfft,Fs);
%     CC1(:,k) = temp;
% end

% amf_ex=reshape(CC1(find(freq==f),:),16,3);

amf_ex=reshape(data_freq(64001,:),16,3);

figure
plot(abs(data_freq(:,1)).')

%% Simulation
% k=1;
% for Theta=linspace(0,2*pi,100)
%     for kk=1:15
%         cost(kk,k)=sum(sum(abs(angle(Gw*exp(-i*w))-angle(circshift(amf_ex*exp(i*Theta),kk)).')));
%     end
%     k=k+1;
% end
% 
% figure
% surf(cost)
% 
% find(cost==min(min(cost)))
% 
% cost(15,780/15)
% 
% Theta(52)

Theta=3.2368;
t=0
figure;
for time=1:length(t)
    subplot(2,1,1)
    imagesc([z1],theta,angle(Gw*exp(-i*w*t(time))));
    title("理论")
    subplot(2,1,2)
    imagesc([z1],theta,angle(circshift(amf_ex*exp(i*Theta),14)).');
    title("试验")

    pause(0.01) 
end
