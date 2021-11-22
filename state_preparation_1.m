% just to understand how much clock time we need I am adding a tstart and
% tend it is not all required but just for the fun.
tstart=tic;
% Also to clear all the stored variables I am giving a clear all but this
% might slow the code a bit but not concerned at the moment.
clear all
% This code basically takes the initial streched state |9/2> to any state
% which we can choose
% j is the spin which is 9/2 for the example at hand.
% J=N( the dimension of the Hilbert space)
% The control Hamiltonian is H(t)=\Omega\left[\cos(\phi(t))I_x+\sin(\phi(t))I_y\right]+\beta I_z^2
j=9/2;
J=2*j+1;
% The algorithm of the code the function Joperators create the I_x, I_y,I_z
% and X_initial and X_final are the initial and final states  and we are
% considering the cat state to be the final state in the case

X_initial=zeros(10,1);
X_initial(10)=1;

X_final=zeros(10,1);
X_final(2)=1/sqrt(2);
X_final(9)=1/sqrt(2);

% time is the total time required for the control to work and in this code
% we are taking time multiplied by \pi so t=3 means time=3\pi
time=3;
beta=1;
Omega=1;


[Jx,Jy,Jz,Jminus,Jplus]= Joperators(j); % creates the angular momentum operators


Lower_limt=-ones(J,1); % basically used for the constrained optimization fmincon

Initial_seed=zeros(J,1); % The initial seed used for the optimization


% The following option allow us to test whether our gradient is accurate as
% well as to consider scenario in which we don't need a gradient or not.

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

% PWC_phase is the function in which we have written both the gradient and
% the function to minimize.

[QQ,fval] = fmincon(@(x)PWC_Phase(x,X_initial,X_final,j,time*pi,beta,Omega),Initial_seed,[],[],[],[],Lower_limt,-Lower_limt,[],options);

% The following steps helps in writting both the final state created and
% control waveform to a single text file 
[Fidelity,Gradient,phase,dt,psi_f]=PWC_Phase(QQ,X_initial,X_final,j,time*pi,beta,Omega);


psi_f=psi_f(:,:,end);
phase=transpose(phase);
Gradient;
dt;
phase1=Initial_seed;
Fidelity;

H=[real(psi_f);imag(psi_f);real(X_final);imag(X_final);phase;Fidelity];
% Writes everything into a single data for convenience.

%writematrix(H,"D:\Graphics_state_prep\sina_"+p+"_beta_"+beta+"_number_"+ii+"_time_"+time+"_control_"+M+".txt",'Delimiter','tab');










