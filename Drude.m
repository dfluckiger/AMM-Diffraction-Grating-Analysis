%       nk = Drude(name,lambda)
%       return vector (real and complex part) for Index of Refraction using
%       Drude model.
%       Input
%	   name:	text, name of material (see below)
%	  lambda: free space wavelength (microns)
%
%	   Copyright (c) 2020 David Fluckiger
%       Author: David Fluckiger
%       Date: May 2020, revision: 1.0
%
%       This code is provided to encourage the development of AMM method.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to 
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 

function nk = Drude(name,lambda)
% given wavelength (lambda) in microns, return index of refraction for
% Drude model of the following material
% ==== NOTE =====
%  lambda is unnormalized, and has units of microns
switch lower(name)
	case 'aluminum'
		p1 = 647;
		p2 = 119000;
	case 'copper' 
		p1 = 278;
		p2 = 63800;
	case 'gold'
		p1 = 216;
		p2 = 72500;
	case 'lead'
		p1= 1450;
		p2 = 62000;
	case 'silver'
		p1 = 145;
		p2 = 72500;
	case 'tungsten'
		p1 = 433;
		p2 = 48400;
	otherwise
		p1 = -100;
		p2 = -100;
end

lambda = 1e4/lambda;		% lambda has units of microns (unnormalized!)
if p1>0
	e1 = p2^2/(lambda^2 + p1^2);
	e2 = p2^2*p1/(lambda*(lambda^2 + p1^2));
	nk = sqrt(e1+1i*e2);
	nk = [real(nk),imag(nk)];
else
	nk = [1,0];
end