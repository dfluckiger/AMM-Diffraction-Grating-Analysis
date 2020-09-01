%       plotEigEq
%       generates a (log(abs(f(mu))+1) contour plot of the eigenvalue equation
%
%	   Copyright (c) 2020 David Fluckiger
%       Author: David Fluckiger
%       Date: May 2020, revision: 1.0
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to 
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
%	Usage
% plotEigEq(ul,lr,Gparms)
%	Input:
%	ul, lr: complex numbers the give upper-left, and lower-right corners
%	of the bounding plot region
%	Gparams:	grating layer defintion structure 
%		Gparms: structure array that defines various parameters specific
%		to a grating layer
%	     Gparams.k0  (=2*pi/lambda)
%		Gparams.kx0  (=k0*sind(theta))
%		Gparams.TMflag  (=true for TM mode, false for TE mode)
%		Gparams.di  (list of region widths)
%		Gparams.mpN (number of decimal digits for MP calculations)
function plotEigEq(ul,lr,Gparams)
% plot EigEq
% set limits on mu

mur = linspace(real(ul),real(lr),100);
mui = linspace(imag(ul),imag(lr),500);
% function
F = zeros(length(mur),length(mui));

for k1 = 1:length(mur)
    for k2 = 1:length(mui)
        mu = mur(k1)+1i*mui(k2);
        [fY,~] = EigEq(mu,Gparams);
        F(k1,k2) = log(abs(fY)+1);
    end
end

contour(mui,mur,F,50);
ax = gca;
ax.YDir = 'normal';

