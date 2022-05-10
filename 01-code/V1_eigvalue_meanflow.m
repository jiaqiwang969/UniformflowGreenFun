clc
clear
close all

%% Add Subfunction
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));


%% Mode Generator
r0=0.5;
w=5;
m = [-20:20];              
n = [1:6];
M=0.5; 
[Base] = BaseJ1(m,n(end)); 
beta=sqrt(1-M^2);
Eigp=(-w*M+sqrt(w^2-beta^2*Base.jmn_pm.^2))/beta^2;  % right running
Eigm=(-w*M-sqrt(w^2-beta^2*Base.jmn_pm.^2))/beta^2;  % left running


Omagp=w-Eigp*M;
Omagm=w-Eigm*M;

Qp = +(Eigp + Omagp*M).*(1-(m./Base.jmn_pm).^2); %R-AIAA-17
Qm = -(Eigm + Omagm*M).*(1-(m./Base.jmn_pm).^2); %R-AIAA-17



for km=1:length(m)
    pp_miu = besselj(m(km),Base.jmn_pm*r0).*besselj(m,alpha*r)./alpha./Q./besselj(m,alpha).^2.*exp(-1i*kapa*xo);
end
