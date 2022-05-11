function [Gnm,Tnm]=greenfun_bf_wall(Tm,Gm,Nxh,i_r,circle_number,mode_enlarge,k_nm,r_pole)   
%% Ϊ�����bf�������ֹ�ʽ���ٶȺ���,135s
Tmz=[];Gmz=[];  

for nk=1:length(mode_enlarge)
    for kk=1:circle_number
        Gmz{nk,kk}= Gm{nk,1}*exp(sqrt(-1)*mode_enlarge(nk)*((kk-1)*2*pi/circle_number));%��������ǶԵĹ�����ͬ�Ƕ�ͨ�����ְ취����Ҫ�ٴ����¼���
        Tmz{nk,kk}= Tm{nk,1}*exp(sqrt(-1)*mode_enlarge(nk)*((kk-1)*2*pi/circle_number));% %�ֱ��������������k�Բ���ֵ
    end   
end

%% 
for kk=1:circle_number
       T{kk,1}=repmat(chebfun(0,[0,2*pi]),1,length(r_pole));
       G{kk,1}=repmat(chebfun(0,[0,2*pi]),1,length(r_pole));
       
    for kkk=1:length(mode_enlarge) 
       G{kk,1}=G{kk}+Gmz{kkk,kk};%��Ϊ�����٣��ۼ����в���
       T{kk,1}=T{kk}+Tmz{kkk,kk};%��Ϊ�����٣��ۼ����в���
   
    end
    G1{kk,1}=reshape(G{kk,1}(linspace(0,2*pi*(1-1/Nxh),Nxh)),Nxh,i_r)';
    T1{kk,1}=reshape(T{kk,1}(linspace(0,2*pi*(1-1/Nxh),Nxh)),Nxh,i_r)';
end
Gnm=cell2mat(G1);
Tnm=cell2mat(T1);

end
