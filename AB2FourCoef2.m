%       AB2FourCoef2
%       Given a set of AB matrices (boundary condition solutions for the
%       layer internal boundaries), compute the Fourier Coefficients for
%       the eigenfunctions
%       
% Usage:
%	[an,kan] = AB2FourCoef2(N,xi,eigCoe,kx0,mpN)
%		N:  Fourier index: -N:N for a total 2*N+1 coefficients
%		xi: vector of index transistions
%		eigCoe: structure array that defines various intermediate values
%		for eigenvalue computation (obtained from call to EigEq)
%		kx0: k0*n1*sind(theta), x-component of incident k-vector
%		mpN: if MP is needed,  number of digits (default = 32);
% Returns
%		an:  array of Fourier coefficients for eigenfunction
%		kan: array of Fourier coefficients for permitivity weighted eigenfunction (TM)
%
%	   Copyright (c) 2020 David Fluckiger
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

function [an,kan] = AB2FourCoef2(N,xi,eigCoe,kx0)
% compute -N:N fourier coefficients for E(x) = aonb*Q(X) + Y(x)
% AB is computed from ABY2x2 for a given eigenvalue (beta)

% Kflag is for TM mode K matrix which is Intgral[ u(x)_n exp(i kx(m) x)/eps(x) dx , {x,0,1}]
% TM mode needs to run twice, once with Kflag, false, for eigVec, and then true
% for Kvec.

mpN = eigCoe.mpN;

if ismp(eigCoe.r)
    an = mp(zeros(2*N+1,1),mpN);
else
    an = zeros(2*N+1,1);
end
kan = an;
tau = exp(1i*kx0);
AB = eigCoe.AB;     %need to accumulate AB
for i = 2:length(xi)-1
      AB(:,:,i) = AB(:,:,i-1)*AB(:,:,i);
end
aonb = AB(2,1,end)/(tau-AB(1,1,end));
ei =  eigCoe.ei;

for k = -N:N
    aq = 0;
    ay = 0;
    kaq = 0;
    kay = 0;
    T = [1,0;0,1];
    for j = 2:length(xi)
        r = eigCoe.beta(j-1);
% assumes kxn = kx0 - 2pi(-M:M)
% integral( un(x) exp(-1i*(kx0-2*pi*(-M:m)))
         t = kx0 - 2*pi*k;
        c = eigCoe.c(j-1);
        s = eigCoe.s(j-1);
        d = xi(j)-xi(j-1); 
        
        if abs(t-r) < 1e-9  % check for sin(x)/x condition
            e12 = exp(-1i*r*(xi(j-1)+2*xi(j)));
            e1  = exp(2i*r*xi(j-1));
            e2  = exp(2i*r*xi(j));
            A = T(1,1);
            B = T(1,2);   
            h = e12*(-e1*(B-1i*A*r)+e2*(B-2i*B*d*r+A*r*(-1i+2*r*d)))/(4*r^2);
            aq = aq + h;
            kaq = kaq + h/ei(j-1);
            A = T(2,1);
            B = T(2,2);  
            h = e12*(-e1*(B-1i*A*r)+e2*(B-2i*B*d*r+A*r*(-1i+2*r*d)))/(4*r^2);
            ay = ay + h;
            kay = kay + h/ei(j-1);
        elseif abs(t+r) < 1e-9
            e1 = exp(-1i*r*xi(j-1));
            e1s = exp(2i*r*xi(j-1));
            e2s = exp(2i*r*xi(j));
            A = T(1,1);
            B = T(1,2);   
            h = e1*(-A*r*(1i*e2s+e1s*(-1i-2*r*d))+B*(-e2s+e1s*(1+2i*r*d)))/(4*r^2);
            aq = aq + h;
            kaq = kaq + h/ei(j-1);
            A = T(2,1);
            B = T(2,2);
            h = e1*(-A*r*(1i*e2s+e1s*(-1i-2*r*d))+B*(-e2s+e1s*(1+2i*r*d)))/(4*r^2);
            ay = ay + h;
            kay = kay + h/ei(j-1);
        else
            e1 = exp(1i*t*xi(j-1));
            e2  = exp(1i*t*xi(j));
            e12  = exp(-1i*t*(xi(j-1)+xi(j)));
            A = T(1,1);
            B = T(1,2);   
            h =  e12*(e2*(B+1i*A*t)*r+e1*(-(B+1i*A*t)*r*c+(-1i*B*t+A*r*r)*s))/(r*(r-t)*(r+t));
            aq = aq + h;
            kaq = kaq + h/ei(j-1);
            A = T(2,1);
            B = T(2,2);   
            h = e12*(e2*(B+1i*A*t)*r+e1*(-(B+1i*A*t)*r*c+(-1i*B*t+A*r*r)*s))/(r*(r-t)*(r+t));
            ay = ay + h;
            kay = kay + h/ei(j-1);
        end
        T = AB(:,:,j-1);
    end
    an(k+N+1) = aq*aonb+ay;   
    kan(k+N+1) = kaq*aonb+kay;   
end
an = double(an);
kan = double(kan);
