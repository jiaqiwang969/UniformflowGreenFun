% TodoList:
% 1. change the normal unit from nomalized unit
% 2. add soft boundary
% 3. uniform flow


format long
format compact
clc;
clear all
close all

%% parameters
r0 = 0.5;          %rotating radius
omegar = 1*pi; %1.1*pi;   %rotating angular speed
omega0 = 1;

rr= linspace(0.001,1,10);
thetar = linspace(0,2*pi,50);

% observer location
%
% r = 0.8;
%theta = pi/2;
for k1=1:length(rr)
    for k2=1:length(thetar)

        xo = 0;
        r = rr(k1);%0.8;
        theta = thetar(k2);%pi/2;

        %% generate the retarded time steps
        T = 2*pi/omegar;
        N = 100;
        dtau = T/N;
        tau = [0:dtau:(10*N-1)*dtau];

        %% observer time step, iteration over m, n
        M=0.0;
        t = zeros(1,10*N);
        mm = [-30:30];
        n = [1:1:10];

        NT = length(t);
        NM = length(mm);
        NMIU = length(n);

        Xs_1 = 0*zeros(1,NT);
        Xs_2 = r0*cos(omegar*tau);    %source position at tau
        Xs_3 = r0*sin(omegar*tau);

        Xo_1 = xo;
        Xo_2 = r*cos(theta);          %observer position at t
        Xo_3 = r*sin(theta);

        t=tau+sqrt((Xo_1-Xs_1).^2 + (Xo_2-Xs_2).^2 + (Xo_3-Xs_3).^2);    % the observer time;

        %% Eigbalue_meanflow.m
        for m1=1:NM

            m = mm(m1);
            omegam = omega0 + m*omegar;
            [Base] = BaseJ1(m,n(end));
            beta=sqrt(1-M^2);
            Eig=((-omegam*M+sign(xo)*sqrt(omegam^2-beta^2*Base.jmn_pm.^2))/beta^2); % R-Fundament-83
            Omag=omegam-Eig*M;

            Q = abs(besselj(m,Base.jmn_pm)).^2.*(1-(m./Base.jmn_pm).^2);

            %     sign(xo) * ((Eig+Omag*M).*(1-m^2./Base.jmn_pm.^2));


            p_n = besselj(m,Base.jmn_pm*r0).*besselj(m,Base.jmn_pm*r)./Base.jmn_pm./Q./besselj(m,Base.jmn_pm).^2.*exp(-1i*Eig*xo);
            p_m = sum(p_n)*exp(-1i*m*theta);
            pp(m1,:) = 1i*exp(1i*(omega0+m*omegar)*t)* p_m;
        end

        p(:,k1,k2)=1/2/pi*sum(pp,1);



    end
end

[Theta,Rho]=meshgrid(thetar,rr);
[yy,xx]=pol2cart(Theta,Rho);




figure
for time=1:1000
    offset=0.05; cont=22; % contour setting
    s1=subplot(1,2,1); contour(xx,yy,reshape(real(p(time,:,:)),size(xx)),cont); ...
        axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
    s2=subplot(1,2,2); contour(xx,yy,reshape(imag(p(time,:,:)),size(xx)),cont); ...
        axis('square'); xlabel(''); ylabel('imag', 'FontSize', 20);
    % sgtitle(['Duct Mode-m',num2str(m),'-n',num2str(n)], 'FontSize', 30)
    pause(0.01) %in seconds

end















