
function [Gv_sum,Tv_sum]=greenfun_bf_ref(Tv,Gv,x0,zH_model,Nxh,i_r,i_theta,k_nm,r_pole)   
%Gnm,Tnm,
%为不同的zH_model赋值
Tvz=[];Gvz=[];
for kz=1:length(zH_model)    
     for nk=1:length(k_nm)%可以并向   
%             Tmz{kz}{nk,1}=Tm{nk,1}*exp(sqrt(-1)*k_nm(nk)*(zH_model(kz)-x0));
%             Gmz{kz}{nk,1}=Gm{nk,1}*exp(sqrt(-1)*k_nm(nk)*(zH_model(kz)-x0)); 
         for kr=1:i_r
            Tvz{kz}{kr,nk}=Tv{1,nk}{kr,1}*exp(sqrt(-1)*k_nm(nk)*(zH_model(kz)-x0));
            Gvz{kz}{kr,nk}=Gv{1,nk}{kr,1}*exp(sqrt(-1)*k_nm(nk)*(zH_model(kz)-x0));
         end
     end
end


%% 为正向的显示
for kz=1:length(zH_model)  %40
 for nkk=1:i_r %30
         temp1{kz}{nkk,1}=zeros(size(Gvz{kz}{1, 1}));
         temp2{kz}{nkk,1}=zeros(size(Tvz{kz}{1, 1}));

      for nk=1:length(mode_enlarge)
         temp1{kz}{nkk,1}=temp1{kz}{nkk,1}+Gvz{kz}{nkk,nk};  
         temp2{kz}{nkk,1}=temp2{kz}{nkk,1}+Tvz{kz}{nkk,nk};  

      end
 end
 
Gv_sum{kz}=temp1{kz};  %30个不同r值点下的，单极子声源,格林函数
Tv_sum{kz}=temp2{kz};  %30个不同r值点下的，偶极子声源,格林函数
end

end