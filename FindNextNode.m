function [ Index ] = FindNextNode(NodesCoord,PrevNode,RefNode,TempNodes)
%  FindNextNode: finds the next node in the candidate region boudary
%                process. The next one (after the reference one) is picked
%                from the fixed set of nodes.
%
% INPUTS
%
%  NodesCoord     : nodes coordinates
%  PrevNode       : previous node
%  RefNode        : reference (current) node
%  TempNodes      : set of nodes
%
% OUTPUTS
%
%  Index          : index of the next node
%
% Copyright (c) 2018 Gdansk University of Technology
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to 
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 
% Author: Piotr Kowalczyk
% Project homepage: https://github.com/PioKow/GRPF
%
P=NodesCoord(PrevNode,:);
S=NodesCoord(RefNode,:);
N=NodesCoord(TempNodes,:);


NoOfTempNodes=size(N,1);


SP=ones(NoOfTempNodes,1)*(P-S);
%SN=N-S

SN=N-ones(size(N,1),1)*S;

LenSP=sqrt(SP(:,1).^2+SP(:,2).^2);
LenSN=sqrt(SN(:,1).^2+SN(:,2).^2);


DotProd=SP(:,1).*SN(:,1)+SP(:,2).*SN(:,2);


Phi=acos(DotProd./(LenSP.*LenSN));

Temp=find(SP(:,1).*SN(:,2)-SP(:,2).*SN(:,1)<0);

Phi(Temp)=2*pi-Phi(Temp);

[~,Index]=min(Phi);



end
