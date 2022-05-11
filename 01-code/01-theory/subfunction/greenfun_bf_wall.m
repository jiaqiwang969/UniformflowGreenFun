function [Gnm,Tnm]=greenfun_bf_wall(Tm,Gm,Nxh,i_r,circle_number,mode_enlarge,k_nm,r_pole)   
%% 为逆向的bf构建格林公式，速度很慢,135s
Tmz=[];Gmz=[];  

for nk=1:length(mode_enlarge)
    for kk=1:circle_number
        Gmz{nk,kk}= Gm{nk,1}*exp(sqrt(-1)*mode_enlarge(nk)*((kk-1)*2*pi/circle_number));%这个代码是对的哈，不同角度通过这种办法不需要再次重新计算
        Tmz{nk,kk}= Tm{nk,1}*exp(sqrt(-1)*mode_enlarge(nk)*((kk-1)*2*pi/circle_number));% %分别求得所以轴向波数k对测点的值
    end   
end

%% 
for kk=1:circle_number
       T{kk,1}=repmat(chebfun(0,[0,2*pi]),1,length(r_pole));
       G{kk,1}=repmat(chebfun(0,[0,2*pi]),1,length(r_pole));
       
    for kkk=1:length(mode_enlarge) 
       G{kk,1}=G{kk}+Gmz{kkk,kk};%仍为无量纲，累加所有波数
       T{kk,1}=T{kk}+Tmz{kkk,kk};%仍为无量纲，累加所有波数
   
    end
    G1{kk,1}=reshape(G{kk,1}(linspace(0,2*pi*(1-1/Nxh),Nxh)),Nxh,i_r)';
    T1{kk,1}=reshape(T{kk,1}(linspace(0,2*pi*(1-1/Nxh),Nxh)),Nxh,i_r)';
end
Gnm=cell2mat(G1);
Tnm=cell2mat(T1);

end
