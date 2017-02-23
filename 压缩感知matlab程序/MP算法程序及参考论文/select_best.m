function [proj,scale,translation,freq,phase]=select_best(signal_r,N,a_base,j_min,j_max,u_base,p_min,v_base,k_min,w_base,i_min,i_max); 

% this subroutine is to select in the dictionary the best atom suited the siganl or the residual of the signal 

% INPUT 
% the signal_r: the signal or the residual of the signal to be decomposed 
% the N: the longth of the signal or of the residual of the signal or the length of the atoms 

% 输入
% the signal_r: 剩余的待分解信号
% N：信号，剩余信号，原子的长度

% parameters :the parameter to construct the dictionary , it have much influence on the speed of the decomposition 
% a_base=2 
% j_min=0; 
% j_max=log2(N); 
 
%OUTPUT 
% proj: the projection of the signal or the residual of the signal on the best atom 
% the scale: the scale of the best atom (s in the formula) 
% the translation : the translation of the best atom (u in the formula) 
% the freq: the frequency of the best atom (v in the formula) 
% phase: the phase of the best atom (w in the formula) 

% 输出
% proj: 在最佳原子基上的剩余信号
% the scale: 最佳原子基的伸缩 
% the translation : 最佳原子基的平移
% the freq: 最佳原子基的频率
% phase: 最佳原子基的相位

% proj_trans :to determine which projection is biggest 
% proj_trans ：判断哪个结果是最好的

proj_trans=0; 
proj=0; 
  
for j=j_min:j_max 
   for p=p_min:N*2^(-j+1) 
      for k=k_min:2^(j+1) 
         for i=i_min:i_max 
            s=a_base^j; % 计算伸缩的值，2为基底 
            u=p*s*u_base; % 计算平移的值
            v=k*(1/s)*v_base; % 计算频率值
            w=i*w_base; % 计算相位值
            t=0:N-1;
            t=(t-u)/s; % 进行伸缩平移变换
             
            g=(1/sqrt(s))*exp(-pi*t.*t).*cos(v*t+w); % 对单个原子信号进行伸缩平移和调制
             
            g=g/sqrt(sum(g.*g)); % 作归一化
             
            proj_trans=sum(signal_r.*g); % 生成重建信号
             
            if abs(proj_trans)>abs(proj) % 判断重建信号是否是最佳（如果是则将相应参数值赋给相应变量）
               proj=proj_trans; 
               scale=s; 
               translation=u; 
               freq=v; 
               phase=w; 
            end 
         end 
      end 
   end 
end 