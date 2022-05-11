%eigvalueInitial-2010-1-5
function [Tm,Tv,Gm,Gv]=greenfun_bf_initial(N,D,r,Ratio,Mx,w,Nxh,i_r,i_theta,mode_enlarge,EigValue,x_pole,r_pole,Angle) %w给定
Entropy=0;Boundary=[1]; Type={'Hard Wall';'Lined Outer Wall';'Lined Inner Wall';'Lined Outer&Inner Wall'};
z_t=1-2*sqrt(-1);z_h=1-2*sqrt(-1);beta=[0.3]; 
Tr=0.0;Omag=0.0;
Tr=0;Omag=0;M_theta=Tr./r+Omag*r; %当旋流系数过大，则必须计算其影响;Mx不能为0
[c02,rou0,P0,s0]=entropyPara(r,N,Ratio,Omag,Tr,Entropy,beta);  

%% 输入参数

for nk=1:length(mode_enlarge)%可以并向
%     [V,lam]=eigfun_AB(r,D,N,w,mode_enlarge(nk),Ratio,Mx,M_theta,rou0,P0,c02,Boundary,z_t,z_h);%求特征
% %     figure;plot(real(lam),imag(lam),'.');hold on;set(gca, 'XLim',[-400 350]); set(gca, 'YLim',[-500 500]);hold on;set(gca,'Color',[1,1,1]); title(['m=',num2str(m(nk)),'  w=',num2str(w)]);
%     crLayer=[min((w-mode_enlarge(nk)*M_theta./r)./Mx);max((w-mode_enlarge(nk)*M_theta./r)./Mx)];
%     cutOffLine=GMM_Cluster3(lam,crLayer); 
%     [mode1,len{1,nk}]=eig_choose_nmax(V,N,r,Ratio,lam,w,Mx,Tr,Omag,-1,10100,-1,500,crLayer,cutOffLine,n_len,0); %修改过bug
%     disp(['m=',num2str(mode_enlarge(nk)),'  w=',num2str(w),'  upstram   mode_enlarge:  ',num2str(len{1,nk})]);
    %区别：Tm,Gm只包含外圈的infomation（用于逆向波数形成），Gm和Gv包含完全的全界面的数据，用于正向仿真显示
    [Tm{nk,1},Tv{1,nk},Gm{nk,1},Gv{1,nk}]=greenfun_dipole_initial(Boundary,mode_enlarge(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,EigValue(nk),z_t,z_h,r_pole,x_pole,i_r,i_theta,Angle,0);%算green函数（1*length(x_pole) cell）
    %[Tm2{nk,1}]=greenfun_dipole(Boundary,mode_enlarge(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,lam(mode1),z_t,z_h,r_pole,x_pole,0,0);%算green函数（1*length(x_pole) cell）

end


end