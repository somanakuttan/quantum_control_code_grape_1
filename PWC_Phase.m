function [Fidelity,Gradient,phase,dt,psi_f] = PWC_Phase(const_phases,psi,psit,spin,t_final,beta1,Omega)
% Piecwise Constant Phase: This function is called by PI_search.m to determine the fidelity and
% gradient with control parameters as phases.
% wf is the number of step sizes

% read in operators and control fields
%dt = const.dt;
dim = 2*spin + 1;%const.dim;
n_iso = 1;%const.n_iso;
% make operators
[Jx,Jy,Jz] = make_fs(spin);
steps_length = length(const_phases);
dt = t_final/steps_length;

% initializes unitaries to be used in GRAPE algorithm
U = zeros(dim,dim,steps_length);
U(:,:,1) = eye(dim);
dU_dphi = zeros(dim,dim,steps_length);
psi_f = zeros(dim,n_iso,(steps_length + 1));
psi_f(:,:,1) = psi;
phase = zeros(1,steps_length);

% iterate over each time step
for jj = 1:steps_length
    phase(jj) = const_phases(jj);
    H_tot = beta1*(Jz^2) + Omega*(cos(pi*phase(jj))*Jx + sin(pi*phase(jj))*Jy);
    % H_tot = H_tot - trace(H_tot)*eye(dim)/dim;
    
    % eigendecomposition of H_tot
    [VH,DH] = eig(H_tot);
    D_vec = diag(DH);
    VH_ct = ctranspose(VH);
    
    exp_eig_factor = -1i*dt*diag(exp(-1i*dt*D_vec));
    
    % solve for derivative of U w.r.t control fields (in H_tot eigenbasis)
    for kk = 1:dim
        for ll = (kk + 1):dim
            exp_eig_factor(kk,ll) = ( exp(-(1i)*D_vec(kk)*dt) - exp(-(1i)*D_vec(ll)*dt) ) / ( D_vec(kk)-D_vec(ll) );
            exp_eig_factor(ll,kk) = exp_eig_factor(kk,ll);
        end
    end
    exp_eig_factor(isnan(exp_eig_factor))=1;
    
    % derivative of H with respect to each control field
    dH_dphi = pi*Omega*(-Jx*sin(pi*phase(jj))+ Jy*cos(pi*phase(jj)));
    
    dU_dphi(:,:,jj) = VH * ((VH_ct*dH_dphi*VH).*exp_eig_factor) * VH_ct;
    
    U(:,:,jj) = (VH*diag(exp(-1i*dt*D_vec))*VH_ct);
    psi_f(:,:,jj + 1) = U(:,:,jj)*psi_f(:,:,jj);
end

trace_term = conj(psit'*psi_f(:,:,steps_length + 1));
Fidelity = - abs(trace_term)^2/n_iso^2;

% solve for gradient
if nargout > 1
    Gradient = zeros(steps_length,1);
    psi_p = zeros(n_iso, dim, steps_length + 1);
    psi_p(:,:,steps_length + 1) = psit';
    
    % solve for psi_p's
    for jj = 1:steps_length
        jp = steps_length - jj + 1;
        psi_p(:,:,jp) = psi_p(:,:,jp + 1)*U(:,:,jp);
        
        % create control wave form update method
        Gradient(jp,1) = - (2/n_iso^2)*real(psi_p(:,:,jp+1)*dU_dphi(:,:,jp)*psi_f(:,:,jp)*trace_term);
    end
end

