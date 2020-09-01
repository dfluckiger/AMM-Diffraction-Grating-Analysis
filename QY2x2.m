%       QY2x2
%       Compute the piecewise-eigenfunctions (Theta, Psi: Q, Y) and then
%       the field over some vector of x-values, good for diagnostic, not
%       used for DE evaluation
%       
% Usage:
%	[Q,Y,Qp,Yp,E,ex] = QY2x2(x,xi,kx0,eigCoe,mpN)
%		x:  vector of x values (0:...:1)
%		xi: vector of index transistions
%		kx0: k0*n1*sind(theta), x-component of incident k-vector
%		eigCoe: structure array that defines various intermediate values
%		for eigenvalue computation (obtained from call to EigEq)
%		mpN: if MP is needed,  number of digits (default = 32);
% Returns
%		Q:  Theta(x)
%		Y:   Psi(x)
%		Qp:  Theta'(x)
%		Yp: Psi'(x)
%		E:  E(x) (for TE) or H(x) (for TM)
%		ex: ei(x) is a piecewise function of the permitivity 
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
function [Q,Y,Qp,Yp,E,ex] = QY2x2(x,xi,kx0,eigCoe)

n = length(eigCoe.ei);
mpN = eigCoe.mpN;

en = ones(n,1);
if eigCoe.TMflag
    if n>1
        en = [1,eigCoe.ei(2:end)./eigCoe.ei(1:end-1)];
    end
end

beta = eigCoe.beta;
if ismp(beta)
    AB = mp(zeros(2,2,n),mpN);
else
    AB = zeros(2,2,n);
end
AB(:,:,1) = [1,0;0,1];
Y = [1,0;0,1];
for i = 2:n+1
    r = beta(i-1)*eigCoe.di(i-1);
    S = [cos(r), -beta(i-1)*sin(r);sin(r)/beta(i-1),cos(r)];
    Y(1,2) = Y(1,2)*en(i-1);
    Y(2,2) = Y(2,2)*en(i-1);
    Y = Y*S;
    AB(:,:,i) = Y;
end

es = 1;
if eigCoe.TMflag
    es = en(1)/en(end);
end
EigEq =Y(1,1)+es*Y(2,2)-2*cos(kx0);

if ismp(eigCoe.r)
    Q = mp(zeros(length(x),1),mpN);
else
    Q = zeros(length(x),1);
end
Y = Q;
E = Q;
Qp = Q;
Yp = Y;
ex = ones(length(x),1);  %this is eps(x), piecewise constant
en = ones(length(eigCoe.ei)+1,1);
if eigCoe.TMflag
    if length(eigCoe.ei)>1
        en = [1,eigCoe.ei(2:end)./eigCoe.ei(1:end-1),1];
    end
end

tau = exp(1i*kx0);
aonb = AB(2,1,end)/(tau-AB(1,1,end));

for i = 1:length(x)
    x1 = x(i);
    n = sum(x1>=xi(1:end-1));
    r = beta(n)*(x1-xi(n));
    S = [cos(r), -beta(n)*sin(r); sin(r)/beta(n),cos(r)];
    U = squeeze(AB(:,:,n));
    U(1,2) = U(1,2)*en(n);
    U(2,2) = U(2,2)*en(n);
    U = U*S;
    Q(i) = U(1,1);
    Y(i) = U(2,1);
    Qp(i) = U(1,2);
    Yp(i) = U(2,2);
    E(i) = aonb*Q(i)+Y(i);
    if eigCoe.TMflag
        ex(i) = eigCoe.ei(n);
    end
end

