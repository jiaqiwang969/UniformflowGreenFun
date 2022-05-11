%ҶƬ����-2019-10-16-wjq
%��ƽ������:50%�ҳ���leading edgeΪ��׼,��ʽ�����Ŀ��Ϊ��������ģ����׼��
%��������Ҷ����Ҷ��ǰԵ���������ϵ������stagger��lean��sweep��������,
%reference��NASA/CR-2001-210762,Hanson,Appendix A

% clc
% clear
% close all
%warning('off');
function h=MATLAB4geomTurbo(xr,yr,zr,xh,yh,zh,xr1,yr1)
h=figure;hold on; axis equal;view([180 -30]);

plot3(xr*1000,yr*1000,zr*1000-170,'.');hold on
plot3(xh*1000,yh*1000,zh*1000-170,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',4,...
    'Marker','diamond',...
    'LineStyle','none');axis equal
plot3(xr1*1000,yr1*1000,zr*1000-170,'.');axis equal
IGV_stagger=0;IGV_lean=0;IGV_sweep=0;
rotor_stagger=45/180*pi;rotor_lean=-15/180*pi;rotor_sweep=30/180*pi;%���Թ���
stator_stagger=-30/180*pi;stator_lean=0/180*pi;stator_sweep=0/180*pi;
%[fname,location]=uigetfile({'*.geomTurbo';'*.txt';'*.dat';'*.*'},'r');%ѡ���ļ�
fname='DTS-02.geomTurbo';
location='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\GreenFunction4Beforming_With_InitialEigenValue'

Data=ReadgeomTurbo(fname)   %��ȡ�ļ�����ȡ��Ҫ������Ϣ����ȡ����data����dataת��cell��ʽ�����洢�����
Data.IGV_suction_rotate=rotateIGV(Data.IGV_suction,[0,1,0],IGV_stagger,[0,0,-104.015]); % rotation the IGV
Data.IGV_pressure_rotate=rotateIGV(Data.IGV_pressure,[0,1,0],IGV_stagger,[0,0,-104.015]);%%���ֶ���˳ʱ����ת����׼Ϊ0��



grid_index=10;
[Y_IGV,Z_IGV] = meshgrid(linspace(Data.IGV_suction_rotate{1, 1}(1,2),Data.IGV_suction_rotate{1, end}(1,2),grid_index),...
    linspace(Data.IGV_suction_rotate{1,1}(1,3),Data.IGV_suction_rotate{1,1}(end,3),grid_index));
X_IGV = zeros(size(Y_IGV));
plane_IGV=surf(X_IGV,Y_IGV,Z_IGV);
[Y_rotor,Z_rotor] = meshgrid(linspace(Data.rotor_suction{1, 1}(1,2),Data.rotor_suction{1, end}(1,2),grid_index),...
    linspace(Data.rotor_suction{1,10}(1,3)-3,Data.rotor_suction{1,10}(1,3)+60,grid_index));
X_rotor = zeros(size(Y_rotor));
plane_rotor=surf(X_rotor,Y_rotor,Z_rotor);
[Y_stator,Z_stator] = meshgrid(linspace(Data.stator_suction{1, 1}(1,2),Data.stator_suction{1, end}(1,2),grid_index),...
    linspace(Data.stator_suction{1,6}(1,3),Data.stator_suction{1,6}(end,3)+2,grid_index));
X_stator = zeros(size(Y_stator))+Data.stator_suction{1,6}(1,1);
plane_stator=surf(X_stator,Y_stator,Z_stator);
set(plane_IGV,'Visible','on');set(plane_rotor,'Visible','on');set(plane_stator,'Visible','on');

Arc_s1=[828.893086090 452.272727273 0]';Arc_e1=[817.5 420.5 0]';O1=[867.5 420.5 0]';
Arc_s2=[1037.370421204 550.7 0]';Arc_e2=[828.893086090 452.272727273 0]';O2=[1037.370421204 280.7 0]';
L1=[1048.5 550.7 0]';L2=[1037.37 550.7 0]';
X1=arcPlot(Arc_s1,Arc_e1,O1);X2=arcPlot(Arc_s2,Arc_e2,O2);
X=[L1 L2 X2(:,1:end-3) X1(:,2:end-1)].';

%plot3(X(1,:),X(2,:),X(3,:));


Surf_bolt=rotsurf([X(:,3) X(:,2)-min(X(:,2)) X(:,1)-max(X(:,1))+min(Data.hub(:,1))],[0 0 1],[0 0 0],linspace(0,2*pi,37),@surf);
set(Surf_bolt, 'FaceColor',[0.68 0.92 1]);

Wall1=[185 0 -170-1400;185 0 -170];
%Surf_wall1=rotsurf(Wall1,[0 0 1],[0 0 0],linspace(-0.8*pi,0.8*pi,37)-1/2*pi,@surf);
Wall2=[sin(linspace(0+1*pi,2/4*pi+1*pi,20))'*120+305 zeros(20,1) cos(linspace(0+1*pi,2/4*pi+1*pi,20))'*120-1400-170];
% Surf_wall2=rotsurf(Wall2,[0 0 1],[0 0 0],linspace(-0.8*pi,0.8*pi,37)-1/2*pi,@surf);
% set(Surf_wall1, 'FaceColor',[0.7 0.92 1]);
% set(Surf_wall2, 'FaceColor',[0.7 0.92 1]);




%set(Surf_wall1, 'FaceColor',[0.68 0.92 1]);
%set(Surf_wall2, 'FaceColor',[0.68 0.92 1]);

Surf_shround=rotsurf([Data.shround(:,2) zeros(size(Data.shround,1),1) Data.shround(:,1)],[0 0 1],[0 0 0],linspace(-0.8*pi,0.8*pi,37)-1/2*pi,@surf);
Surf_hub=rotsurf([Data.hub(:,2) zeros(size(Data.hub,1),1) Data.hub(:,1)],[0 0 1],[0 0 0],linspace(-1*pi,1*pi,37)-1/2*pi,@surf)
set(Surf_shround, 'FaceColor',[0.93 0.70 0.13]);
set(Surf_hub,'LineStyle','-.',...
    'FaceColor',[0.16 0.38 0.27],...
    'EdgeColor',[0.86 0.86 0.86]);

[X11]=yexing2surf(Data.IGV_suction_rotate);
[X12]=yexing2surf(Data.IGV_pressure_rotate);
IGV(1)=surf(X11{1},X11{2},X11{3},'EdgeColor',[1 0 0]);hold on
IGV(2)=surf(X12{1},X12{2},X12{3},'EdgeColor',[1 0 0]);
Rotor{1}=plotCellData(Data.rotor_suction,'r');hold on;
Rotor{2}=plotCellData(Data.rotor_pressure,'r');hold on;
Stator{1}=plotCellData(Data.stator_suction,'r');hold on;
Stator{2}=plotCellData(Data.stator_pressure,'r');hold on;
set(IGV,'Visible','on');
set(Rotor{1},'Visible','on');set(Rotor{2},'Visible','on')
set(Stator{1},'Visible','on');set(Stator{2},'Visible','on');

%rotate_blade(IGV,[0 0 1],32,'surf'); 
%rotate_blade(Rotor{1},[0 0 1],29,'plot'); rotate_blade(Rotor{2},[0 0 1],29,'plot'); 
%rotate_blade(Stator{1},[0 0 1],37,'plot'); rotate_blade(Stator{2},[0 0 1],37,'plot'); 

%��������ϵ
[para.XAXIS(1) para.XAXIS(2) para.XAXIS(3)]=arrow3d( [0 0 0.5]',3,90,8,'x'); 
[para.YAXIS(1) para.YAXIS(2) para.YAXIS(3)]=arrow3d( [0 0 0.5]',3,90,8,'y');
[para.ZAXIS(1) para.ZAXIS(2) para.ZAXIS(3)]=arrow3d( [0 0 0.5]',3,90,8,'z');

% para.XAXIS(4)=text(110,0,0,'X','Color','r','FontSize',12);
% para.YAXIS(4)=text(0,110,0,'Y','Color','r','FontSize',12);
para.ZAXIS(4)=text(0,10,120,'Z','Color','r','FontSize',12);
%�����������ϵ
origin_IGV=[0 Data.IGV_suction_rotate{1, 1}(1,2:3)];
origin_rotor=[0 Data.rotor_suction{1, 1}(1,2) Data.rotor_suction{1,10}(1,3)-3];
origin_stator=[Data.stator_suction{1, 1}(1,1) Data.stator_suction{1, 1}(1,2) Data.stator_suction{1,6}(1,3)];
% 
% [para.XAXIS_IGV]=arrow3d_relative(origin_IGV',30,8,1.5); 
% [para.XAXIS_rotor]=arrow3d_relative(origin_rotor',30,8,1.5);
% [para.XAXIS_stator]=arrow3d_relative(origin_stator',30,8,1.5);

IGV_new=surf_newxyz(X_IGV,Y_IGV,Z_IGV,origin_IGV,IGV_stagger,IGV_lean,IGV_sweep,grid_index);
rotor_new1=surf_newxyz(X_rotor,Y_rotor,Z_rotor,origin_rotor,rotor_stagger,0,0,grid_index);
rotor_new2=surf_newxyz(X_rotor,Y_rotor,Z_rotor,origin_rotor,rotor_stagger,rotor_lean,0,grid_index);
rotor_new2=surf_newxyz(X_rotor,Y_rotor,Z_rotor,origin_rotor,rotor_stagger,rotor_lean,rotor_sweep,grid_index);
stator_new=surf_newxyz(X_stator,Y_stator,Z_stator,origin_stator,stator_stagger,stator_lean,stator_sweep,grid_index);

function New=surf_newxyz(X,Y,Z,origin,stagger,lean,sweep,grid_index)
%ע�⣺��תΪ�����������ϵ��������ת��Ϊ��������ϵ
newxyz =([X(:), Y(:), Z(:)]-repmat(origin,length(X(:)),1))*Q_sweep(sweep)*Q_lean(lean)*Q_stagger(stagger)+repmat(origin,length(X(:)),1);%ע�������˵��Ⱥ�˳��
New=surf(reshape(newxyz(:,1),grid_index,grid_index),reshape(newxyz(:,2),grid_index,grid_index),reshape(newxyz(:,3),grid_index,grid_index),'EdgeColor',[1 1 1]);
end



function output=Q_stagger(theta)
output=[ cos(theta)    0  -sin(theta) ;
          0            1     0
       sin(theta)     0  cos(theta);];%A-2:��y����ʱ��Ϊ��
end
function output=Q_lean(theta)
output=[  cos(theta)   sin(theta)  0;
          -sin(theta)   cos(theta)   0 ;
            0            0          1];%����stagger�任��A-4:��z����ʱ��Ϊ��   
end
function output=Q_sweep(theta)
output=[  1    0             0 ;
          0  cos(theta)  sin(theta);
          0  -sin(theta)  cos(theta);];%����stagger��lean�任��A-6:��x����ʱ��Ϊ��   
end





%% subfunction
function [X]=yexing2surf(Blade)
for k=1:length(Blade)
s(k)=size(Blade{1, k},1);
end
s1=sort(s);
s_max=s1(end-1);
s_min=min(s);
xx1=[];xx2=[];xx3=[];
for k=1:length(Blade)
    xyz=[];
    xyz=resample(Blade{1,k},s_min,size(Blade{1,k},1));
    xx1=[xx1 xyz(:,1)];
    xx2=[xx2 xyz(:,2)];
    xx3=[xx3 xyz(:,3)];
end
X={xx1;xx2;xx3};
end

function H=plotCellData(h,color)
   [m,n]=size(h);
   for k=1:n
    M=h{k};
    H(k)=plot3(M(:,1),M(:,2),M(:,3),color);hold on
   end
end 

function h1=rotateIGV(h,azel,alpha,origin) % ��תIGV�����մ�ǰԵ��βԵ��˳���������׼Ϊ0��


% find unit vector for axis of rotation
    u = azel(:)/norm(azel);


alph = (alpha+198.2)*pi/180;  %ģ�͵ĽǶ�Ϊ-18.2��
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);y = u(2);z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
   
   [m,n]=size(h); %newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
   
   for k=1:n
    M=h{k};
    [m1,n1]=size(M);
    B=repmat(origin,m1,1);
    M1 = (M-B)*rot+B;
    M2=(fliplr(M1'))';%������ת�ӽ�180�ȣ�ǰβԵ����
    h1{k}=M2;   
   end
  
end

function Data=ReadgeomTurbo(fname)
      ZR_mark=0;%��־hub��shround��
      Z_mark=0;%��־IGV_leading=1;IGV_trailing=2;rotor_leading=3;rotor_trailing=4;stator_leading=5;stator_trailing=6;
      mark=0;%��־IGV_suction-1;IGV_pressure-2;rotor_suction-3;rotor_pressure-4;stator_suction-5;stator_pressure-6;
      fid=fopen(fname,'r');     
      
    for j=1:10000  %��j��
        if ~feof(fid)  %ֻҪ�ĵ�������
            
            line=fgetl(fid);%��һ��һ�ж���ȥ
               if (~isempty(strfind(line,'NI_BEGIN zrcurve')))  %��ȡhub��shround��
                ZR_mark=ZR_mark+1;              
                line=fgetl(fid);
                ZR_numberOfPoints(ZR_mark) = fscanf(fid, '%d',1);  %��ȡĳЩ������Ϣ
                for n = 1:ZR_numberOfPoints(ZR_mark)
                         fscanf(fid, '%c',1);
                         z=fscanf(fid, '%f',1);
                         r=fscanf(fid, '%f',1);   
                         fscanf(fid, '%c',1);    
                         if ZR_mark==1
                         Z_hub(n,1)= [z];
                         R_hub(n,1)= [r];
                         else
                         Z_shround(n,1)= [z];
                         R_shround(n,1)= [r];
                         end
                             
                 end
               
                                    switch ZR_mark
                            case 1
                                Data.hub=[Z_hub,R_hub];
                                clear Z_hub
                                clear R_hub
                            case 2
                                Data.shround=[Z_shround,R_shround];
                                clear Z_shround
                                clear Z_shround
                                    end
                                     
                 end                    
%% ��ȡǰԵ��βԵ����

             
                if ((~isempty(strfind(line,'NI_BEGIN nisolid_angle_at_leading_edge')))|(~isempty(strfind(line,'NI_BEGIN nisolid_angle_at_trailing_edge'))))  %��ȡhub��shround��
                  Z_mark=Z_mark+1;              
                  %line=fgetl(fid);
                  Z_numberOfPoints(Z_mark) = fscanf(fid, '%d',1);  %��ȡĳЩ������Ϣ
                  for n = 1:Z_numberOfPoints(Z_mark)
                         fscanf(fid, '%c',1);
                         x=fscanf(fid, '%f',1);
                         y=fscanf(fid, '%f',1); 
                         z=fscanf(fid, '%f',1); 
                         fscanf(fid, '%c',1);    
                         X1(n,:)= [x,y,z];
                         end
                             
                
            
                                    switch Z_mark
                                                                    case 1
                                Data.IGV_leading=[X1];
                            case 2
                                Data.IGV_trailing=[X1];
                            case 3
                                Data.rotor_leading=[X1]; 
                            case 4
                                Data.rotor_trailing=[X1];
                            case 5
                                Data.stator_leading=[X1];
                            case 6
                                Data.stator_trailing=[X1];
                            end
                                    clear X1
                end 
                
%% ��ȡ Ҷ������    
             if (~isempty(strfind(line,'suction')))|(~isempty(strfind(line,'pressure'))) %�����Ƶ���־��ʱ  |strfind(line,'pressure')
                mark=mark+1;              
                line=fgetl(fid);
                section = fscanf(fid, '%d',1);  %��ȡĳЩ������Ϣ
                continue;
             elseif isempty(strfind(line,['#','   section ']))    % else if isempty(strfind(line,['#','   section ',num2str(k)]))
                 temp=0;  %���״̬
             elseif ~isempty(strfind(line,['#','   section ']))
                        temp=1;
                        k=str2num(line(isstrprop(line,'digit')));  %��ȡ�ı��е�����
                        line=fgetl(fid);
                        numberOfPoints(mark,k) =fscanf(fid, '%d',1);
                        
                       for n = 1:numberOfPoints(mark,k)
                         fscanf(fid, '%c',1);
                         x=fscanf(fid, '%f',1);
                         y=fscanf(fid, '%f',1);
                         z=fscanf(fid, '%f',1);
                         fscanf(fid, '%c',1);
                        
                         X(n,:)= [x,y,z]; 
                       end
                       
                        switch mark
                            case 1
                                Data.IGV_suction{k}=[X];    
                            case 2
                                Data.IGV_pressure{k}=[X];
                            case 3
                                Data.rotor_suction{k}=[X];
                            case 4
                                Data.rotor_pressure{k}=[X];
                            case 5
                                Data.stator_suction{k}=[X];
                            case 6
                                Data.stator_pressure{k}=[X];
                        end
                        clear X  %���ܼӷ��
                        
                 end
             end
       
    end
      Data.numberOfPoints=numberOfPoints; 
    end
                                
function h=rotsurf(curve,direct,point,theta,f)
% rotsurf(curve,dirct,orgin,alpha,fun)���ڻ�����ת����
%   curve=[x,y,z]Ϊĸ�ߣ�����x,y,zΪ���������ֱ����ĸ�ߵ���ά����
%   direct��origin�ֱ������ת��ķ���͸���ת���ϵ�����һ������꣬����������������ȷ����һ��ֱ�ߣ�����ת�ᣬ���У�
%       direct��ʾ��ת��ķ��������ֱ�ʾ��[theta,phi]��[x0,y0,z0]�����У�
%           theta������xoyƽ���x����������ʱ����ת�Ļ��ȣ�phi�����xoyƽ����z����������ת�Ļ���
%           [x0,y0,z0]����������
%           directĬ��Ϊ[0 0 1]����z�᷽��
%       origin=[xo,yo,zo]Ϊ����ת���ϵ�����һ�����꣬Ĭ��Ϊ[0 0 0]��ԭ��
%   ����alphaΪ��ת�Ļ��ȣ�Ĭ��Ϊ0:pi/1:2*pi��������ķ�Χ���ܶȶ������ֶ�����
%   funΪ��ͼ���������Ĭ��Ϊ@mesh
% h=rotsurf(...)
%  ���������ͬʱ���ظ�����ľ��h
%
%��1������ĸ��Ϊx=0,y^2+z^2=1����ת��Ϊx=1,z=-y-2��Բ��
%t=linspace(-pi,pi,37).'; 
%y=sin(t);z=cos(t);x=y-y; 
%rotsurf([x,y,z],[0 -1 1],[1 -2 0],[],@surf) 
%xlabel('x');ylabel('y');zlabel('z');axis equal
%��2������ĸ��Ϊz=x,y=1����ת��Ϊz��ĵ�Ҷ˫����
% t=linspace(-pi,pi,37).'; 
% x=t;y=t-t+1;z=t;
% rotsurf([x y z]) 
% xlabel('x');ylabel('y');zlabel('z');axis equal

assert(nargin>=1 && nargin<=5,'�������������뿴������');
if nargin<5
    f=@mesh;
    if nargin<4
        theta=linspace(0,2*pi,37);
        if nargin<3
            point=[0,0,0];
            if nargin<2
                direct=[0,0,1];
            end
        end
    end
end
curve=squeeze(curve);
assert(ismatrix(curve),'����1��ʽ�����뿴������');
if size(curve,2)~=3
    curve=curve.';
end
assert(size(curve,2)==3,'����1��ʽ�����뿴������');
direct=squeeze(direct);
if isempty(direct)
    direct=[0,0,1];
end
assert(numel(direct)==2 || numel(direct)==3,'����2��ʽ�����뿴������');
if numel(direct)==2
    direct=[cos(direct(2))*cos(direct(1)),cos(direct(2))*sin(direct(1)),sin(direct(2))];
end
direct=direct/norm(direct);
point=squeeze(point);
if isempty(point)
    point=[0 0 0];
end
assert(numel(point)==3,'����3��ʽ�����뿴������');

theta=squeeze(theta);
if isempty(theta)
    theta=linspace(0,2*pi,37);
end
assert(length(theta)==numel(theta) ,'����4��ʽ�����뿴������');

f=squeeze(f);
if isempty(f)
    f=@mesh;
end
if ischar(f)
    assert(length(f)==numel(f),'����5��ʽ�����뿴������');
    f=str2func(f);
end
assert(numel(f)==1 && isa(f,'function_handle'),'����5��ʽ�����뿴������');


x0=point(1);
y0=point(2);
z0=point(3);

x=curve(:,1);
y=curve(:,2);
z=curve(:,3);

nx=direct(1);
ny=direct(2);
nz=direct(3);

[X,~]=meshgrid(x,theta);
[Y,~]=meshgrid(y,theta);
[Z,T]=meshgrid(z,theta);

sint=sin(T);
cost=cos(T);

XX=(X-x0).*(nx^2*(1-cost)+cost)+(Y-y0).*(nx*ny*(1-cost)-nz*sint)+(Z-z0).*(nx*nz*(1-cost)+ny*sint)+x0;
YY=(X-x0).*(ny*nx*(1-cost)+nz*sint)+(Y-y0).*(ny^2*(1-cost)+cost)+(Z-z0).*(ny*nz*(1-cost)-nx*sint)+y0;
ZZ=(X-x0).*(nz*nx*(1-cost)-ny*sint)+(Y-y0).*(nz*ny*(1-cost)+nx*sint)+(Z-z0).*(nz^2*(1-cost)+cost)+z0;

hh=f(XX,YY,ZZ);
if nargout==1
    h=hh;
end
end

function h=rotate_blade(h,azel,N,type)
%h=rotate(h,azel,alpha,origin)
%�ں���rotate�Ļ����ϸĽ���ר��������ҶƬ��ת

% Determine the default origin (center of plot box).
  if ~ishghandle(h)
    error(message('MATLAB:rotate:InvalidHandle'));
  end
  ax = ancestor(h(1),'axes');
  if isempty(ax) || ax==0,
    error(message('MATLAB:rotate:InvalidHandle'));
  end
  matlab.ui.internal.UnsupportedInUifigure(ancestor(ax,'figure'));
  origin = sum([get(ax,'xlim')' get(ax,'ylim')' get(ax,'zlim')'])/2;

% find unit vector for axis of rotation
if numel(azel) == 2 % theta, phi
    theta = pi*azel(1)/180;
    phi = pi*azel(2)/180;
    u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
elseif numel(azel) == 3 % direction vector
    u = azel(:)/norm(azel);
end

for k=1:N-1
alph = 2*pi/N*k;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

for i=1:numel(h),
  t = get(h(i),'type');
  skip = 0;
  if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'patch')
    
    % If patch, rotate vertices  
    if strcmp(t,'patch')
       verts = get(h(i),'Vertices');
       x = verts(:,1); y = verts(:,2); 
       if size(verts,2)>2
          z = verts(:,3);
       else
          z = [];
       end
       
    % If surface or line, rotate {x,y,z}data   
    else
       x = get(h(i),'xdata');
       y = get(h(i),'ydata');
       z = get(h(i),'zdata');
    end
    
    if isempty(z)
       z = -origin(3)*ones(size(y));
    end
    [m,n] = size(z);
    if numel(x) < m*n
      [x,y] = meshgrid(x,y);
    end
  elseif strcmp(t,'text')
    p = get(h(i),'position');
    x = p(1); y = p(2); z = p(3);
  elseif strcmp(t,'image')
    x = get(h(i),'xdata');
    y = get(h(i),'ydata');
    z = zeros(size(x));
  else
    skip = 1;
  end
  
  if ~skip,
    [m,n] = size(x);
    newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
    newxyz = newxyz*rot;
    newx = origin(1) + reshape(newxyz(:,1),m,n);
    newy = origin(2) + reshape(newxyz(:,2),m,n);
    newz = origin(3) + reshape(newxyz(:,3),m,n);
    
    if strcmp(type,'surf') 
        surf(newx,newy,newz,'EdgeColor',[1 1 1]);hold on
        elseif strcmp(type,'plot')
        plot3(newx,newy,newz,'Color',[1 1 1]);hold on
    end
        
        

%     if strcmp(t,'surface') || strcmp(t,'line')
%       surf(newx,newy,newz,'EdgeColor',[1 0 0]);hold on
% %     set(h(i),'xdata',newx,'ydata',newy,'zdata',newz);
%     elseif strcmp(t,'patch')
%       set(h(i),'Vertices',[newx,newy,newz]); 
%     elseif strcmp(t,'text')
%       set(h(i),'position',[newx newy newz])
%     elseif strcmp(t,'image')
%       set(h(i),'xdata',newx,'ydata',newy)
%     end
  end
end
end
end

function [ hHead,hTail,hPlate ] = arrow3d( orgin,r,h1,h2,dir)
%ARROW3D Summary of this function goes here
%  Detailed explanation goes here

% 
% orgin=[0 0 0.2]';
% r=1;
% h=90;
% dir='y';

n=50;
% p=0.15;

hTail=Cylinder(orgin,r/2,h1,dir,n,'open');
hold on

if dir=='x'
    org=orgin+[h1; 0; 0];
elseif dir=='y'
    org=orgin+[0; h1; 0];
elseif dir=='z'
    org=orgin+[0; 0; h1];
end

xlabel('x')
ylabel('y')
zlabel('z')

[hHead hPlate]=Cone( org,r,h2,dir,n,'closed' );
set(hTail,'EdgeAlpha',0)
end

function [ Arrow ] = arrow3d_relative(orgin,h1,r1,r2)
    Arrow(1)=arrow(orgin,orgin+[h1; 0; 0],r1,'BaseAngle',45,'Width',r2,'Color','g');
    Arrow(2)=arrow(orgin,orgin+[0; h1+40; 0],r1,'BaseAngle',45,'Width',r2,'Color','g');
    Arrow(3)=arrow(orgin,orgin+[0; 0; h1],r1,'BaseAngle',45,'Width',r2,'Color','g');
end


function [hCylinder,hPlate1,hPlate2] = Cylinder( orgin,r,h,dir,n,closed )
% This Function plots a cylinder at specified orgin, with specified radius, height, direction(axis), no of points along the circumference
% 
%   Typical Call Cylinder( orgin,r,h,dir,n ), Note: there is a matlab function "cylinder"
% 
%   orgin: vector of order 3 x 1 specifies orgin
%   h    : Height of the Cylinder
%   dir  : String to specify axis of extrution 'x' or 'y' or 'z' (Only along these axes...!)
%   n    : no of points along the circumference
% 
   

t=linspace(0,2*pi,n)';

x1=r*cos(t);
x2=r*sin(t);
h1=orgin(3);
h2=h1+h;


if dir=='y'
    xx1=[[x1;x1(1)] [x1;x1(1)]]+orgin(1);
    xx2=[repmat(h1,length(x1)+1,1) repmat(h2,length(x1)+1,1)]+orgin(2);
    xx3=[[x2;x2(1)] [x2;x2(1)]]+orgin(3);
elseif dir =='x'
    xx1=[repmat(h1,length(x1)+1,1) repmat(h2,length(x1)+1,1)]+orgin(1);
    xx2=[[x1;x1(1)] [x1;x1(1)]]+orgin(2);
    xx3=[[x2;x2(1)] [x2;x2(1)]]+orgin(3);
else
    xx1=[[x1;x1(1)] [x1;x1(1)]]+orgin(1);
    xx2=[[x2;x2(1)] [x2;x2(1)]]+orgin(2);
    xx3=[repmat(h1,length(x1)+1,1) repmat(h2,length(x1)+1,1)]+orgin(3);
end
hCylinder=surf(xx1,xx2,xx3,repmat(3,size(xx1)));

if strcmp(closed,'closed')==1
    hold on
    hPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'g');
    hPlate2=fill3(xx1(:,2),xx2(:,2),xx3(:,2),'g');
end
end

function [ hCone, hPlate1] = Cone( orgin,r,h,dir,n,closed )
%CONE Summary of this function goes here
%  Detailed explanation goes here

x=[r 0]';
y=[0 h]';

[hCone hPlate1 hplate2]=Revolve(orgin,x,y,dir,n,closed);
end

function [hRevolve, hPlate1,hPlate2]=Revolve(orgin,x,y,dir,n,closed)

m=length(x);
t=linspace(0,2*pi,n);
if dir=='y'
    xx1=x*sin(t)+orgin(1);
	xx2=y*ones(1,n)+orgin(2);
	xx3=x*cos(t)+orgin(3);    
elseif dir =='x'
    xx1=y*ones(1,n)+orgin(1);
	xx2=x*sin(t)+orgin(2);
	xx3=x*cos(t)+orgin(3);
else
    xx1=x*cos(t)+orgin(1);
	xx2=x*sin(t)+orgin(2);
	xx3=y*ones(1,n)+orgin(3);
    
end
hRevolve=surf(xx1,xx2,xx3,repmat(2,size(xx1)));

if strcmp(closed,'closed')==1
    hold on
    hPlate1=fill3(xx1(1,:),xx2(1,:),xx3(1,:),'g');
    hPlate2=fill3(xx1(m,:),xx2(m,:),xx3(m,:),'g');
end
end
%% д���ĵ���
function modify(DATA,IGV_angle)

[fname,location]=uigetfile({'*.geomTurbo';'*.txt';'*.dat';'*.*'},'r');
%copyfile('DTS-02.geomTurbo', 'modify6.geomTurbo');  %�ȿ���һ����Ϊ����,modify��Ϊ�޸��ļ�
    fid=fopen(fname,'r');
    fd=fopen([['New_angle_',num2str(IGV_angle),'_'] ,fname],'w');

%words=input('Give the string you want to delete:\n','s');
%words='47';%��־��

     mark=0;%���mark=1��IGV_suction����mark=2��IGV_pressure����mark=3��rotor_suction����mark=4��rotor_pressure����mark=5��stator_suction����mark=6��stator_pressure����
    
    for j=1:10000  %��j��
        if ~feof(fid)  %ֻҪ�ĵ�������
           
            line=fgetl(fid);%��һ��һ�и�����ȥ
             if (~isempty(strfind(line,'suction')))|(~isempty(strfind(line,'pressure'))) %�����Ƶ���־��ʱ  |strfind(line,'pressure')
                         mark=mark+1;
                fprintf(fd,'%s\r\n',line);  %�����ڸ��Ƹ���Ϣ 
                line=fgetl(fid);
                fprintf(fd,'%s\r\n',line); 
                section = fscanf(fid, '%d',1);  %��ȡĳЩ������Ϣ
                fprintf(fd,'%s',num2str(section));  %�����ڸ��Ƹ���Ϣ    %%���������ṹ���Դ��������
                continue;
             elseif isempty(strfind(line,['#','   section ']))    % else if isempty(strfind(line,['#','   section ',num2str(k)]))
                 temp=0;  %���״̬
                 fprintf(fd,'%s\r\n',line);
             else if ~isempty(strfind(line,['#','   section ']))
                        temp=1;
                        fprintf(fd,'%s\r\n',line);
                        k=str2num(line(isstrprop(line,'digit')));  %��ȡsection number'
                        line=fgetl(fid);
                        fprintf(fd,'%s\r\n',line);
                        numberOfPoints(k) =fscanf(fid, '%d',1);
                        fprintf(fd,'%s\r',num2str(numberOfPoints(k)));   
                        line=fgetl(fid);
                        fprintf(fd,'%s\r\n',line);
                        
                        if mark==1
                             matrix=DATA.IGV_suction_rotate{1,k};  %�滻������
                        elseif mark==2
                             matrix=DATA.IGV_pressure_rotate{1,k};  %�滻������
                        elseif mark==3
                             matrix=DATA.rotor_suction{1,k};  %�滻������
                       elseif mark==4
                             matrix=DATA.rotor_pressure{1,k};  %�滻������
                       elseif mark==5
                             matrix=DATA.stator_suction{1,k};  %�滻������
                       elseif mark==6
                             matrix=DATA.stator_pressure{1,k};  %�滻������
                        end
                        
                        [m,n]=size(matrix);
                            if m~=numberOfPoints(k)
                                disp('����������')
                                break;
                             end
                        for kk=1:numberOfPoints(k)
                            temp=2;
                            line=fgetl(fid);
                            %fprintf(fd,'%s\r\n',line);  %�滻��
                            
                                fprintf(fd,'%f %f %f\n',matrix(kk,:)); 
                            
                        end    
                                
                               
                            
                        end
                            
                      
               %%���������Ե��Ա�֤�����ǻ��޷��뵽��Ч�������ݵĽṹ
                           
                         end
                     
                   
                 
             

             
        else
            break;
        end
    end
     fclose(fd);
    if(temp>10)
        delete(['G' cell2mat(fname(i))]);
    end
    fclose(fid);
    disp('�����ɹ�����')
   
end


function X=arcPlot(Arc_s,Arc_e,O)
R1=Arc_s-O;
R2=Arc_e-O;
theta=subspace(R1,R2);
R=[];
R=[R1,R];
X=[];
W=cross(R1,R2);
th=linspace(0,theta,20);
dt=(th(2)-th(1))/norm(W);
for m=1:length(th)
    Rota=R(:,end)+cross(W,R(:,end))*dt;
    R=[R,Rota];
    X=[X,O+Rota];
end
X=[Arc_s,X,Arc_e];
end

end
