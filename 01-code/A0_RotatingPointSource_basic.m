% this origin code is by ZhaoHuan Wang
% TodoList:
% 1. change the normal unit from nomalized unit
% 2. add soft boundary
% 3. uniform flow 


format long
format compact
clc;
clear all
%% parameters
r0 = 0.5;          %rotating radius
omegar = 1.1*pi;   %rotating angular speed
omega0 = 1;
% observer location
xo = 1;
r = 0.8;
theta = pi/2;

%% generate the retarded time steps
T = 2*pi/omegar;
N = 100;
dtau = T/N;
tau = [0:dtau:(10*N-1)*dtau];

%% observer time step, iteration over m, miu
t = zeros(1,10*N);
mm = [-30];
Miu = [1:1:100];

NT = length(t);
NM = length(mm);
NMIU = length(Miu);

Xs_1 = 0*zeros(1,NT);
Xs_2 = r0*cos(omegar*tau);    %source position at tau
Xs_3 = r0*sin(omegar*tau);

Xo_1 = xo;
Xo_2 = r*cos(theta);          %observer position at t
Xo_3 = r*sin(theta);

t=tau+sqrt((Xo_1-Xs_1).^2 + (Xo_2-Xs_2).^2 + (Xo_3-Xs_3).^2);    % the observer time;

%% Eigbalue_meanflow.m
M=0.0;
n=[1:100]; 


%% calculation 
for m1=1:NM
   
    m = mm(m1);
    omegam = omega0 + m*omegar;


[Base] = BaseJ1(m,n(end)); 
beta=sqrt(1-M^2);
Eigp=(-omegam*M+sqrt(omegam^2-beta^2*Base.jmn_pm.^2))/beta^2;  % right running
Eigm=(-omegam*M-sqrt(omegam^2-beta^2*Base.jmn_pm.^2))/beta^2;  % left running
Omagp=omegam-Eigp*M;
Omagm=omegam-Eigm*M;
Qp = +(Eigp + Omagp*M).*(1-(m./Base.jmn_pm).^2); %R-AIAA-17
Qm = -(Eigm + Omagm*M).*(1-(m./Base.jmn_pm).^2); %R-AIAA-17

%     alpha0 = 1/2*(Miu+1/2*m-3/4)*pi+1/2*sqrt((Miu+1/2*m-3/4).^2*pi^2-2*m^2+1/2);
%     alpha = fsolve(@(alpha)(besselj(m-1,alpha)-besselj(m+1,alpha))/2,alpha0,optimoptions('fsolve','Display','off'));
%     kapa = -i*sqrt(-omegam^2+alpha.^2) %right-running modes
%     
% figure
% plot(kapa,'b.')
% hold on
% plot(real(Eigp),imag(Eigp),'b o')
% plot(real(Eigm),imag(Eigm),'r o')


    Q = kapa./alpha.*(1-(m./alpha).^2);
    pp_miu = besselj(m,alpha*r0).*besselj(m,alpha*r)./alpha./Q./besselj(m,alpha).^2.*exp(-1i*kapa*xo);
    pp_m = sum(pp_miu)*exp(-1i*m*theta);
    
    
pp(m1,:) = 1i*exp(1i*(omega0+m*omegar)*t)* pp_m;
end

p=1/2/pi*sum(pp,1);

figure
h=plot(t,real(p),'-r');
set(h,'LineWidth',0.5);
grid on

xlabel('time (s)');
ylabel('pressure (Pa)');













