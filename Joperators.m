function [Jx,Jy,Jz,Jminus,Jplus]= Joperators(j)
% JOPERATORS Quantum mechanical angular momentum operators.
%
%   [Jx,Jy,Jz,Jminus,Jplus] = JOperators(j)
%
% The operators are with respect to a basis which ranges  from <j -j| to <j
% j| in increments of m=1.  For instance, in the spin-1/2 case, the first
% basis vector is spin-down, and the second is spin-up. This is somewhat
% contrary to Physics convention.
% 
% Factors of Planck's constant have been omitted from all operators, and
% the raising and lowering operators (J+ and J-) have an implicit square
% root.  This is to keep the results exactly representable.  To correct for
% this, something like the following is necessary:
%
% Jx = hbar * Jx;                  and same for Jy and Jz
% Jminus = hbar * Jminux .^ 0.5;   and same for Jplus
%
% Tobin Fricke <tobin@splorg.org>
% University of Rochester 2005-01-30

mvalues = -j:1:j;
cardinality = length(mvalues);

% Jminus is zero except for directly above the diagonal
% There is an implicit \hbar and square root

Jminus = zeros(cardinality);
for i=2:cardinality,
    m = mvalues(i);
    Jminus(i-1, i) = sqrt(j*(j+1) - m*(m-1));
end

Jplus = zeros(cardinality);
for i=1:(cardinality-1),
    m = mvalues(i);
    Jplus(i+1, i) = sqrt(j*(j+1) - m*(m+1));
end

i = sqrt(-1);

Jx = (1/2) * (Jplus + Jminus);
Jy = (i/2) * (Jplus - Jminus);
Jz = diag(mvalues);
end

