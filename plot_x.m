clear all
beta=1;
Omega=1;
time=7;

% A=importdata("D:\Graphics_state_prep\cat_state_11.txt");
A=importdata("D:\Graphics_state_prep\xgate_17.txt");
 QQ = A(1:110)/pi;
%QQ=A(21:50);
X0=zeros(10,1);
X0(10)=1;
X1=zeros(10,1);
X1(1)=1;
AA=eye(10);
for i=1:10
    X0=AA(:,i);
[Fidelity,phase,dt,psi_f]=PWC_Phase_extension(QQ,X0,X1,9/2,time*pi,beta,Omega,5);
        Fidelity
        
        
    
        
     
        
  writematrix(real(psi_f),"D:\Graphics_state_prep\final_x_real_"+i+".txt",'Delimiter','tab');
    writematrix(imag(psi_f),"D:\Graphics_state_prep\final_x_imag_"+i+".txt",'Delimiter','tab');
end
u = repelem(QQ,5);
writematrix(u,"D:\Graphics_state_prep\final_x_wave.txt",'Delimiter','tab')
