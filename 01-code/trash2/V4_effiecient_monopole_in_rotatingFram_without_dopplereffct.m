format long
format compact
clc;
clear all
tic
%% parameters
r0 = 0.1; %rotating radius
omegar = 50*pi; %rotating angular speed
omega0 = 10;
% observer location
        rr= linspace(0.001,1,10);
        thetar = linspace(0,2*pi,50);

for k1=1:length(rr)
    for k2=1:length(thetar)
        xo = 0.5;
        r = rr(k1);%0.8;
        theta = thetar(k2);%pi/2;

        %% generate the retarded time steps
        T = 2*pi/omegar;
        N = 100;
        dtau = T/N;
        tau = [0:dtau:(10*N-1)*dtau];

        %% observer time step, iteration over m, miu
        t = zeros(1,10*N);
        mm = [-30:30];
        Miu = [1:1:10];

        NT = length(t);
        NM = length(mm);
        NMIU = length(Miu);

        Xs_1 = 0*zeros(1,NT);
        Xs_2 = r0*cos(omegar*tau);   %source position at tau
        Xs_3 = r0*sin(omegar*tau);

        Xo_1 = xo;
        Xo_2 = r*cos(theta);   %observer position at t
        Xo_3 = r*sin(theta);

        t=tau+sqrt((Xo_1-Xs_1).^2 + (Xo_2-Xs_2).^2 + (Xo_3-Xs_3).^2);    % the observer time;

        %% calculation
        for m1=1:NM

            m = mm(m1);
            omegam = omega0 + m*omegar;

            alpha0 = 1/2*(Miu+1/2*m-3/4)*pi+1/2*sqrt((Miu+1/2*m-3/4).^2*pi^2-2*m^2+1/2);
            alpha = fsolve(@(alpha)(besselj(m-1,alpha)-besselj(m+1,alpha))/2,alpha0,optimoptions('fsolve','Display','off'));
            kapa = sqrt(omegam^2-alpha.^2); %right-running modes

            Q = kapa.*(1-(m./alpha).^2);
            pp_miu = besselj(m,alpha*r0).*besselj(m,alpha*r)./alpha./Q./besselj(m,alpha).^2.*exp(-1i*kapa*xo);

            pp_m = sum(pp_miu)*exp(-1i*m*theta);


            pp(m1,:) = 1i*exp(1i*(omega0+m*omegar)*t)* pp_m;
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
toc







