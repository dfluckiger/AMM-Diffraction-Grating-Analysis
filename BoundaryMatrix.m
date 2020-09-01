function [BC,d,mu] = BoundaryMatrix(BC,j,J1,J2,K1,K2,Roots,k1z,k3z,n1,n3,L,Gparams,eFlag,d2,mu2)
%       Boundary 
%       Create Boundary Condition Matrix for multi-layer structure
%
%       Author: David Fluckiger
%       Date: May 2020, revision: 1.0
%
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to 
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 
%       This code is provided to encourage the development of AMM method.
%
% Usage: BC = Boundary(BC,j,J,K,k1z,k3z,n1,n3,N,eFlag))
%
%  outputs:
%         BC    updated BC matrix for current layer
%
% inputs:
%    BC:   matrix to be update of size (N*number_of_layers)^2
%    j:             number of current layer
%    J1,J2:      current and prior layer Fourier-Rayleigh coefficients of eigenfunction of current layer
%    K1,K2:    like J but for TM mode (1/eps(x) weighting in integral)
%    Roots:    eigenvalues/roots
%    k1z:  superstrate kz vector components
%    k3z:  substrate kz vector components
%    n1:   superstrate index of refraction
%    n3:   substrate index of refraction
%    L:      layer thickness
%    Gparams:   grating layer parameters
%    eFlag:    true for last layer, false otherwise
% Notes: 
%

N = size(J1,2);
M = size(J1,1);
k0 = Gparams.k0;
d = diag(exp(1i*k0*Roots*L));
mu = diag(1i*Roots);
if Gparams.TMflag
     ik1 = diag(k0*n1^2./(1i*k1z));
     ik3 = diag(k0*n3^2./(1i*k3z));     
else
     ik1 = diag(k0./(1i*k1z));
     ik3 = diag(k0./(1i*k3z));
     K1 = J1;         % for TE K not used
     K2 = J2;
end

id = 1:N;
ifd = 1:M;
if j==1             %first layer, interface with superstrate
     BC(ifd,id) = J1+ik1*K1*mu;
     BC(ifd,id+N) = J1*d-ik1*K1*mu*d;
end
if eFlag   % last layer, interface with substrate
     ro = M*(2*j-1) + (1:M);
     co = N*(2*j-2) + (1:N);
     BC(ro,co) = J1*d-ik3*K1*mu*d;
     BC(ro,co+N) = J1+ik3*K1*mu;                
end
if j>1
     ro = (2*j-3)*M + (1:M);              % note J1 = current, J2 = previous
     co = (2*(j-2))*N + (1:N);
     BC(ro,co) = J2*d2;
     BC(ro,co+N) = J2;
     
     BC(ro+M,co) = -K2*mu2*d2;
     BC(ro+M,co+N) = K2*mu2;
     
     BC(ro,co+2*N) = -J1;
     BC(ro,co+3*N) = -J1*d;
     
     BC(ro+M,co+2*N) =  K1*mu;
     BC(ro+M,co+3*N) = -K1*mu*d;
end
     
     
