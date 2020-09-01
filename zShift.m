%       zShift
%       Enforce root search to stay inside some box
%       
% Usage:
%	z = zShift(z,Uul,Llr)
%		z:  complex number 
%		Uul, Llr, are complex numbers that define a search volume (prevent walk-off)
%		
% Returns
%		z inside box
%		if root moves beyond boundary, z is the point closest to the
%		boundary
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
function z = zShift(z,Uul,Llr)
% move z to boundary if it is outside box
a = real(z);
b = imag(z);
if a>real(Uul), a = real(Uul); end
if a<real(Llr), a = real(Llr); end
if b>imag(Llr), b = imag(Llr); end
if b<imag(Uul), b = imag(Uul); end
z = a + 1i*b;