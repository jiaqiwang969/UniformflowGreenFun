clc
clear
close all

%% Add Subfunction
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));


%% parameters
f=6000;
w = f*2*pi/343*0.3;
r0 = 0.5; theta0=0; z0=0;  % source term
r = linspace(0,1,200); theta = 0; z1 = linspace(-1,0,200);z1(end)=[];z2 = linspace(0,1,200);  % observer location
M = 0;
Angle=0;
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
    Gmn1(:,:,:,km)=(-Eigm_mn(:,km)*cos(Angle)+m(km)./r0*sin(Angle)).*besselj(m(km),Base.jmn_pm(:,km)*r0)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qm_mn(:,km))...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigm_mn(:,km)*z1)),n,1,length(z1));
    Gmn2(:,:,:,km)=(-Eigp_mn(:,km)*cos(Angle)+m(km)./r0*sin(Angle)).*besselj(m(km),Base.jmn_pm(:,km)*r0)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qp_mn(:,km))...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigp_mn(:,km)*z2)),n,1,length(z2));
end
Gm1=reshape(sum(Gmn1,1),length(r),length(z1),length(m))./(2);       %AIAA-20
Gm2=reshape(sum(Gmn2,1),length(r),length(z2),length(m))./(2);       %AIAA-20
% Gw=Gm*exp(-1i*m*theta).';    %AIAA-20

Gw1=sum(Gm1,3);
Gw2=sum(Gm2,3);
Gw=[Gw1 Gw2];

figure;
t=1:0.01:10;
for time=1:length(t)
    subplot(2,1,1)
    imagesc([z1 z2],r,real(Gw*exp(-i*w*t(time))));
    axis xy;
    axis equal
    subplot(2,1,2)
    imagesc([z1 z2],r,imag ...
        (Gw*exp(-i*w*t(time))));
    axis xy;
    axis equal
    pause(0.01) 
end

