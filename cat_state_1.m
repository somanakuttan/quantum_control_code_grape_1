%This program basically tries to simulate the control waveform for the
%case of the random initial seed and with the initial state being the
%streched state along the z and considered 20 final state
%(20 of them and take the best fidelity out of it).
clear all % clears all previous variables for safety

beta=1;

Omega=1;


j=9/2; % spin under consideration

J=2*j+1; % J is the dimension of the Hilbert space

[Jx,Jy,Jz,Jminus,Jplus]= Joperators(j); % creates the angular momentum operators
M=3*J;
[A,B]=eig(Jz); % creates eigenstates along Jz
X0=A(:,J);
C11=-ones(M,1); % giving the constraint such that the output control
%wave form lies between -1 and 1

    time=2
    
 X1=1/sqrt(2)*(A(:,J-1)+A(:,2)); 
    
    
    
        
        for jj=1:20
            % parallel loop to find the best out of the 20 random initial
            % guesses
            a = -1;
            b = 1;
            C0=(b-a).*rand(M,1) + a;
            % The above three steps create random number between a and b
            
            options = optimset(...  % these are the settings from unitary control search
                'TolX',            1e-16,...    % related to minimum tolerance for change in x
                'TolFun',          1e-8,...    % minumum tolerance for change in the function
                'MaxIter',         2000,...   % maximum number of iterations
                'DerivativeCheck', 'off',...    % compare analytic gradient to numerical estimation (off)
                'GradObj',         'on',...         % tells matlab whether the gradient is supplied to fminunc
                'LargeScale',       'off', ... % when turned off will increase calculation speed
                'Display',         'off',...          % output type in matlab command window, USE 'iter' to turn on
                'MaxFunEvals',     10^6,...         % maximum number of function calls'MaxFunEvals', 500);
                'ObjectiveLimit',  -0.9999);
            
            % options above creates option for optimization
            [QQ,fval] = fmincon(@(x)PWC_Phase(x,X0,X1,j,time*pi,beta,Omega),C0,[],[],[],[],C11,-C11,[],options)
            % creates final contrl wave form and final fidelity
            
            
            F(jj)=-fval;
            F1(:,jj)=QQ;
            F2(:,jj)=C0;
        end
        [M,I] = max(F);
        QQ=F1(:,I);
        [Fidelity,Gradient,phase,dt,psi_f]=PWC_Phase(QQ,X0,X1,j,time*pi,beta,Omega);
        Fidelity
        C0=F2(:,I);
        
        % The above three steps create best out of the 20 initial control
        % seeds and find the fidelity, control wave form and final state.
        
        psi_f1=psi_f(:,:,end);
        phase=transpose(phase);
        Gradient;
        dt;
        phase1=C0;
        
        H=[real(psi_f1);imag(psi_f1);phase;Fidelity;real(X1);imag(X1);phase1];
        % Writes everything into a single data for convenience.
        
        writematrix(H,"D:\Graphics_state_prep\cat_state_11.txt",'Delimiter','tab');
        
  writematrix(real(psi_f),"D:\Graphics_state_prep\final_cat_state_real_1.txt",'Delimiter','tab');
    writematrix(imag(psi_f),"D:\Graphics_state_prep\final_cat_state_imag_1.txt",'Delimiter','tab');

