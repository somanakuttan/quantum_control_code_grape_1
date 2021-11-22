AA=importdata("D:\Graphics_state_prep\xgate_11.txt")/pi;
QQ=AA(1:110);
time=6;
beta=1;
Omega=1;
A=eye(10)
for i=1:10
i
X0=A(:,i);
a=mod(i+1,10)
X1=A(:,3);

[Fidelity,Gradient,phase,dt,psi_f]=PWC_Phase(QQ,X0,X1,9/2,time*pi,beta,Omega);
        Fidelity
        
        
    
        
     
        
  writematrix(real(psi_f),"D:\Graphics_state_prep\final_x_state_real_"+i+".txt",'Delimiter','tab');
    writematrix(imag(psi_f),"D:\Graphics_state_prep\final_x_state_imag_"+i+".txt",'Delimiter','tab');
end
