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
Omega=1;
r0 = 0.8; theta0=0; z0=0;  % source term
r = linspace(0,1,50); theta = linspace(0,2*pi,90); z1 = linspace(-1,0,3);z1(end)=[];z2 = linspace(0,1,3);  % observer location
M = 0;
%% Mode Generator
m = [5];
n = [1];

[Base] = BaseJ1(m,n);
beta=sqrt(1-M^2);

% Omagam_mn = w - M*Eigm_mn;   %AIAA-9
% Omagap_mn = w - M*Eigp_mn;   %AIAA-9
% Qm_mn =  (Eigm_mn+Omagam_mn*M).*(1-m.^2./Base.jmn_pm.^2); %Rienstra
% Qp_mn =  (Eigm_mn+Omagam_mn*M).*(1-m.^2./Base.jmn_pm.^2);

t=0:0.00000001:0.000001;
for km=1:length(m)

   wr =w + m(km)*Omega;
   kappa_mn=sqrt(wr^2-beta^2*Base.jmn_pm(:,km).^2);
   Eigm_mn=(-wr*M+kappa_mn)/beta^2;  % left running
   Eigp_mn=(-wr*M-kappa_mn)/beta^2;  % right running
   Qm_mn =  kappa_mn.*(1-m(km).^2./Base.jmn_pm(:,km).^2); %Lowis
   Qp_mn =  kappa_mn.*(1-m(km).^2./Base.jmn_pm(:,km).^2);

    Gmn1(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r0)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qm_mn)...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigm_mn*z1)),n,1,length(z1));
    Gmn1t(:,:,:,km,:)=bsxfun(@times,Gmn1(:,:,:,km),reshape(exp(-i*(wr)*t),1,1,1,length(t)));
    Gmn2(:,:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*r0)...
        ./besselj(m(km),Base.jmn_pm(:,km)).^2./(Qp_mn)...
        .*besselj(m(km),Base.jmn_pm(:,km)*r)...
        .*reshape(exp(-i*(Eigp_mn*z2)),n,1,length(z2));
    Gmn2t(:,:,:,km,:)=bsxfun(@times,Gmn2(:,:,:,km),reshape(exp(-i*(wr)*t),1,1,1,length(t)));
end

Gm1t=reshape(sum(Gmn1t,1),length(r),length(z1),length(m),length(t))./(-2*pi*i);       %AIAA-20
Gm2t=reshape(sum(Gmn2t,1),length(r),length(z2),length(m),length(t))./(-2*pi*i);       %AIAA-20
tem1=exp(-1i*m.'*(theta-theta0));
tem2=exp(-1i*m.'*(theta-theta0));

Gw1t=reshape(sum(bsxfun(@times,reshape(Gm1t,length(r),length(z1),length(t),length(m)),reshape(tem1,1,1,1,length(m),length(theta))),4),length(r),length(z1),length(theta),length(t));
Gw2t=reshape(sum(bsxfun(@times,reshape(Gm2t,length(r),length(z2),length(t),length(m)),reshape(tem2,1,1,1,length(m),length(theta))),4),length(r),length(z2),length(theta),length(t));

Gw=[Gw1t Gw2t];


[Theta,Rho]=meshgrid(theta,r);
[yy,xx]=pol2cart(Theta,Rho);

figure
for time=1:length(t)
    offset=0.05; cont=22; % contour setting
    s1=subplot(2,2,1); contour(xx,yy,reshape(real(Gw(:,length(z1)+1,:,time)),size(xx)),cont); ...
    axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
    s2=subplot(2,2,2); contour(xx,yy,reshape(imag(Gw(:,length(z1)+1,:,time)),size(xx)),cont); ...
    axis('square'); xlabel(''); ylabel('imag', 'FontSize', 20);
    subplot(2,2,3);imagesc([z1 z2],r,real(Gw(:,:,1,time)));
    axis('square'); axis xy;axis equal;xlabel(''); ylabel('real', 'FontSize', 20);ylim([0 1]);
    subplot(2,2,4);imagesc([z1 z2],r,imag(Gw(:,:,1,time))); 
    axis('square'); axis xy;axis equal;xlabel(''); ylabel('imag', 'FontSize', 20);ylim([0 1]);
    pause(0.01) %in seconds
end




