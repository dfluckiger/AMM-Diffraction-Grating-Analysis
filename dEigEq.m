%       dEigEq
%       Evaluate the derivative eigenvalue equation for the layer defined
%       by previous call to EigEq 
%       
% Usage:
%	dfY = dEigEq(eigCoe)
%		eigCoe is a structure array which holds all the intermediate
%		results from a call to EigEq(mu,Gparams)
% Returns
%		dfY derivative of eigenfunction at mu
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

function dY = dEigEq(eigCoe)

% comput 1st derivative of Y(mu) function (trace of product of AB's)
%  d EigEq/ dmu
% for each section expand derivative using product rule
k0 = eigCoe.k0;
k2 = k0^2;
T =  [0,0;0,0];

for i = 1:eigCoe.n
    dm = eigCoe.di(i)*eigCoe.mu;
    dd = eigCoe.di(i)^2*k2*sqrt(eigCoe.ei(i));
    if abs(eigCoe.beta(i)) == 0.
        S = [dd,2*dd/eigCoe.di(i)*eigCoe.en(i);dd*eigCoe.di(i)/3,dd*eigCoe.en(i)];
    else
        S = [k2*dm*eigCoe.sb(i), k2*(dm*eigCoe.c(i)+eigCoe.mu*eigCoe.sb(i))*eigCoe.en(i);...
               k2*(-dm*eigCoe.c(i)+eigCoe.mu*eigCoe.sb(i))/eigCoe.beta(i)^2, k2*dm*eigCoe.sb(i)*eigCoe.en(i)];
    end

    Y = [1,0;0,1];
     for j = 1:eigCoe.n
        if i==j
            Y = Y*S;
        else
            Y = Y*squeeze(eigCoe.AB(:,:,j));
        end
    end
    T = T+Y;
end

dY = T(1,1)+ eigCoe.es*T(2,2);
