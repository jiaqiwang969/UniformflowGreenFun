clc
clear
close all

%% Add Subfunction
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));
%



%% parameters
f=6000;
w = f*2*pi/343*0.3;
r0= 0.8; theta0=[0:15]/16*2*pi; z0=0; phi0=[0 pi 0 pi 0 pi 0 pi 0 pi 0 pi 0 pi 0 pi];

r = linspace(0,1,20); theta = linspace(0,2*pi,180); z2 =[10];  % observer location
M = 0;
%% Mode Generator
m = [-50:50];
n = [50];

[Base] = BaseJ1(m,n);

beta=sqrt(1-M^2);
kappa_mn=sqrt(w^2-beta^2*Base.jmn_pm.^2);

% Eigm_mn=(-w*M+kappa_mn)/beta^2;  % left running
Eigp_mn=(-w*M-kappa_mn)/beta^2;  % right running


% Omagam_mn = w - M*Eigm_mn;   %AIAA-9
Omagap_mn = w - M*Eigp_mn;   %AIAA-9
% Qm_mn =  kappa_mn.*(1-m.^2./Base.jmn_pm.^2); %Lowis
Qp_mn =  kappa_mn.*(1-m.^2./Base.jmn_pm.^2);

% Qm_mn =  (Eigm_mn+Omagam_mn*M).*(1-m.^2./Base.jmn_pm.^2); %Rienstra
% Qp_mn =  (Eigm_mn+Omagam_mn*M).*(1-m.^2./Base.jmn_pm.^2);


for km=1:length(m)
%     Gmn1(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r0)...
%         ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qm_mn(:,km))...
%         .*besselj(m(km),Base.jmn_pm(:,km)*r)...
%         .*reshape(exp(-i*(Eigm_mn(:,km)*z1)),n,1,length(z1));
    Gmn2(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r0)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qp_mn(:,km))...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigp_mn(:,km)*z2)),n,1,length(z2));
end
% Gm1=reshape(sum(Gmn1,1),length(r),length(z1),length(m))./(-2*pi*i);       %AIAA-20
Gm2=reshape(sum(Gmn2,1),length(r),length(z2),length(m))./(-2*pi*i);       %AIAA-20


t=0:0.01:10;
for mt=1:length(theta0)
   tem0=exp(-i*m.'*(theta-theta0(mt)));
   % Gw01=reshape(sum(bsxfun(@times,Gm1,reshape(tem0,1,1,length(m),length(theta))),3),length(r),length(z1),length(theta));
   Gw02=reshape(sum(bsxfun(@times,Gm2,reshape(tem0,1,1,length(m),length(theta))),3),length(r),length(z2),length(theta));
   tim0=exp(-i*w*t+phi0(mt));
   s02(:,:,:,:,mt)=bsxfun(@times,Gw02,reshape(tim0,1,1,1,length(t)));

end

signal=sum(s02,5);





[Theta,Rho]=meshgrid(theta,r);
[yy,xx]=pol2cart(Theta,Rho);
figure
for time=1:length(t)
    offset=0.05; cont=22; % contour setting
    s1=subplot(2,2,1); contour(xx,yy,reshape(real(signal(:,end,:,time)),size(xx)),cont); ...
    axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
    s2=subplot(2,2,2); contour(xx,yy,reshape(imag(signal(:,end,:,time)),size(xx)),cont); ...
    axis('square'); xlabel(''); ylabel('imag', 'FontSize', 20);
    subplot(2,2,3);imagesc([z2],r,real(signal(:,:,1,time)));
    axis('square'); axis xy;axis equal;xlabel(''); ylabel('real', 'FontSize', 20);ylim([0 1]);
    subplot(2,2,4);imagesc([z2],r,imag(signal(:,:,1,time))); 
    axis('square'); axis xy;axis equal;xlabel(''); ylabel('imag', 'FontSize', 20);ylim([0 1]);
    pause(0.01) %in seconds
end


