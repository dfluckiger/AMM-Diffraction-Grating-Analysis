%  rect_dom: generates the initial mesh for rectangular domain z=x+jy x\in[xb,xe] , y\in[yb,ye]
%
% INPUTS
%
%  xb     : real part range begin 
%  xe     : real part range end 
%  yb     : imag part range begin 
%  ye     : imag part range end 
%  r      : initial mesh step
%
% OUTPUTS
%
%  NewNodesCoord     : generated nodes coordinates
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
function NewNodesCoord = rectdom( xb,xe,yb,ye,r )

nc = max(2,ceil(abs(xe-xb)/r))+1;
nr = max(2,ceil(abs(ye-yb)/r))+1;
X = linspace(xb,xe,nc)';
Y1 = linspace(yb,ye,nr)';
Y2 = Y1 + (Y1(2)-Y1(1))/2;
Y2 = [yb;Y2(1:end-1);ye];
NewNodesCoord = [];
for j = 1:nc
	if mod(j,2) == 1
		xy = [repmat(X(j),[nr,1]),Y1];
		NewNodesCoord = [NewNodesCoord;xy];
	else
		xy = [repmat(X(j),[nr+1,1]),Y2];
		NewNodesCoord = [NewNodesCoord;xy];
	end
end


