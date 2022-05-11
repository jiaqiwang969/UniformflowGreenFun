clc
clear
close all

%% Add Subfunction
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));


%% Mode Generator
m = [0:20];              
n = [1];
[Base] = BaseJ1(m,n(end)); 
M=0.5; 
beta=sqrt(1-M^2);

w=linspace(0,20,1000);
for k=1:length(w)
    Eigp(k,:)=(-w(k)*M+sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2;  % right running
    Eigm(k,:)=(-w(k)*M-sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2;  % left running
end



figure
for k=1:length(m)
   scatter(w,real(Eigp(:,k)),'.'); 
    hold on 
    scatter(w,real(Eigm(:,k)),'.');
end
ylabel({'Re(k_{m1})'});
xlabel({'w'});
title({'Figure-14-fundamental of Duct acoustic% right running'});
ylim([-20,20])
