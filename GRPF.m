function [nrts,rts,rtsMul] = GRPF(NewNodesCoord, Tol, ItMax, NodesMax, SkinnyTriangle, Gparams,bc)
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

% origianl code can be found at the above location, this script has been
% slightly modified for the AMM problem

NodesCoord=[];
NrOfNodes=0;
Quadrants = [];

fun = @(z) EigEq(z,Gparams);
%fun = @(z) (z-1)*(z-1i).^2*(z+1)^3./(z+1i);

r = bc(end)/2;
xb = bc(1)+r; xe = bc(2)-r; yb = bc(3)+r; ye = bc(4)-r; 
Uul = xe + yb*1i;
Llr = xb + ye*1i;
%% general loop
warning('off','all');
it=0;
while it<ItMax && NrOfNodes<NodesMax
    it=it+1;
    NodesCoord=[NodesCoord ; NewNodesCoord];
    
    zFunctionValues = complex(zeros(size(NewNodesCoord,1),1));
    zQuadrants = zeros(size(NewNodesCoord,1),1);
    for Node=1:size(NewNodesCoord,1)
	   x = NodesCoord(NrOfNodes+Node,1);
	   y = NodesCoord(NrOfNodes+Node,2);
        z=x+1i*y;
	   zFunctionValues(Node) = fun(z);
        zQuadrants(Node) = vinq( zFunctionValues(Node,1) );
    end
    Quadrants = [Quadrants;zQuadrants];
    
    NrOfNodes=size(NodesCoord,1);  
    DT = delaunayTriangulation(NodesCoord(:,1),NodesCoord(:,2)); 
    Elements = DT.ConnectivityList;
    Edges = edges(DT); 
    NrOfElements=size(Elements,1);  
    PhasesDiff=mod(Quadrants(Edges(:,1))-Quadrants(Edges(:,2)),4);
    CandidateEdges=Edges(PhasesDiff==2|isnan(PhasesDiff),:);
	
%    vis( NodesCoord, Edges, Quadrants, PhasesDiff)
%    figure(4)
%	triplot(DT)
%	hold on;
%	xy =  [NodesCoord(CandidateEdges(:,1),:),NodesCoord(CandidateEdges(:,2),:)];
%	for i = 1:size(xy,1)
%		px = [xy(i,1);xy(i,3)];
%		py = [xy(i,2);xy(i,4)];
%		line(px,py,'Linewidth',3);
%	end
%
%	for i = 1:length(Quadrants)
%		switch abs(Quadrants(i))
%			case 1
%				plot(NodesCoord(i,1),NodesCoord(i,2),'*r');
%			case 2
%				plot(NodesCoord(i,1),NodesCoord(i,2),'*y');
%			case 3
%				plot(NodesCoord(i,1),NodesCoord(i,2),'*g');
%			case 4
%				plot(NodesCoord(i,1),NodesCoord(i,2),'*b');
%		end
%	end
%	hold off;
%	drawnow;
 
	if isempty(CandidateEdges)
	    rts = [];
	    nrts = 0;
	    rtsMul = [];
        return
    end
    
    Nodes1OfCandidateEdges=CandidateEdges(:,1);
    Nodes2OfCandidateEdges=CandidateEdges(:,2);
    
    CoordinatesOfNodes1OfCandidateEdges=NodesCoord(Nodes1OfCandidateEdges,:);
    CoordinatesOfNodes2OfCandidateEdges=NodesCoord(Nodes2OfCandidateEdges,:);
    
    CandidateEdgesLengths=sqrt(sum((CoordinatesOfNodes2OfCandidateEdges-CoordinatesOfNodes1OfCandidateEdges).^2,2));
    MinCandidateEdgesLengths=min(CandidateEdgesLengths);
    MaxCandidateEdgesLengths=max(CandidateEdgesLengths);
   
    if MaxCandidateEdgesLengths<Tol
	    % should never come here, if so something is wrong (probaly tol is
	    % to small).
        break
    end       

    Temp=CandidateEdgesLengths>Tol;
    ReduCandidateEdges=CandidateEdges(Temp,:);

    Temp=zeros(NrOfNodes,1);
    Temp(ReduCandidateEdges(:,1))=1;
    Temp(ReduCandidateEdges(:,2))=1;
    CandidateNodes=find(Temp==1);
    
    ArrayOfCandidateElements = vertexAttachments(DT,CandidateNodes);
    Temp=zeros(NrOfElements,1);
    for k = 1:size(ArrayOfCandidateElements,1)
        Temp(ArrayOfCandidateElements{k})=Temp(ArrayOfCandidateElements{k})+1;
    end 
    
    IDOfFirstZoneCandidateElements=find(Temp>1);
    IDOfSecondZoneCandidateElements=find(Temp==1);
  
    NoOfFirsrZoneCandidateElements=size(IDOfFirstZoneCandidateElements,1);
    FirstZoneCandidateElements=Elements(IDOfFirstZoneCandidateElements,:);    
         
    for k = 1:NoOfFirsrZoneCandidateElements      
        TempExtraEdges((k-1)*3+1,:)=[FirstZoneCandidateElements(k,1) FirstZoneCandidateElements(k,2)];
        TempExtraEdges((k-1)*3+2,:)=[FirstZoneCandidateElements(k,2) FirstZoneCandidateElements(k,3)];
        TempExtraEdges((k-1)*3+3,:)=[FirstZoneCandidateElements(k,3) FirstZoneCandidateElements(k,1)];      
    end
    
    NewNodesCoord=sum(NodesCoord(TempExtraEdges(1,:),:))/2;
    
    for k = 2:3*NoOfFirsrZoneCandidateElements
        CoordOfTempEdgeNode1=NodesCoord(TempExtraEdges(k,1),:);
        CoordOfTempEdgeNode2=NodesCoord(TempExtraEdges(k,2),:);
        TempNodeCoord=(CoordOfTempEdgeNode1+CoordOfTempEdgeNode2)/2;
        TempEdgeLength=sqrt(sum((CoordOfTempEdgeNode2-CoordOfTempEdgeNode1).^2));
        
        if TempEdgeLength>Tol
            DistNodes=sqrt((NewNodesCoord(:,1)-TempNodeCoord(1)).^2 +(NewNodesCoord(:,2)-TempNodeCoord(2)).^2 );
            if sum(DistNodes<2*eps)==0
                NewNodesCoord=[NewNodesCoord ; TempNodeCoord];
            end
	   end
    end
     
    % removing the first new node if the edge is too short
    CoordOfTempEdgeNode1=NodesCoord(TempExtraEdges(1,1),:);
    CoordOfTempEdgeNode2=NodesCoord(TempExtraEdges(1,2),:);
    TempEdgeLength=sqrt(sum((CoordOfTempEdgeNode2-CoordOfTempEdgeNode1).^2));
    if TempEdgeLength<Tol
        NewNodesCoord(1,:)=[];        
    end
    NoOfSecondZoneCandidateElements=size(IDOfSecondZoneCandidateElements,1);
    SecondZoneCandidateElements=Elements(IDOfSecondZoneCandidateElements,:);         
    for k = 1:NoOfSecondZoneCandidateElements  
        NodesInTempElement=SecondZoneCandidateElements(k,:);
        Node1Coord=NodesCoord(NodesInTempElement(1),:);
        Node2Coord=NodesCoord(NodesInTempElement(2),:);
        Node3Coord=NodesCoord(NodesInTempElement(3),:);
        TempLengths(1)=sqrt(sum((Node2Coord-Node1Coord).^2));
        TempLengths(2)=sqrt(sum((Node3Coord-Node2Coord).^2));
        TempLengths(3)=sqrt(sum((Node1Coord-Node3Coord).^2));
        if max(TempLengths)/min(TempLengths)>SkinnyTriangle
            TempNodeCoord=(Node1Coord+Node2Coord+Node3Coord)/3;
            NewNodesCoord=[NewNodesCoord ; TempNodeCoord];
	   end 
    end
end

% Evaluation of contour edges from all candidates edges 
ArrayOfCandidateElements = edgeAttachments(DT,CandidateEdges(:,1),CandidateEdges(:,2));
Temp=zeros(NrOfElements,1);
for k = 1:size(ArrayOfCandidateElements,1)
    Temp(ArrayOfCandidateElements{k})=1;
end

IDOfCandidateElements=find(Temp==1);
NoOfCandidateElements=size(IDOfCandidateElements,1);
CandidateElements=Elements(IDOfCandidateElements,:);

TempEdges=zeros(NoOfCandidateElements*3,2);

for k = 1:NoOfCandidateElements
    TempEdges((k-1)*3+1,:)=[CandidateElements(k,1) CandidateElements(k,2)];
    TempEdges((k-1)*3+2,:)=[CandidateElements(k,2) CandidateElements(k,3)];
    TempEdges((k-1)*3+3,:)=[CandidateElements(k,3) CandidateElements(k,1)];
end

% Reduction of edges to contour
MultiplicationOfTempEdges=zeros(3*NoOfCandidateElements,1);
RevTempEdges=fliplr(TempEdges);
for k = 1:3*NoOfCandidateElements
    if MultiplicationOfTempEdges(k)==0
        NoOfEdge=find(RevTempEdges(:,1)==TempEdges(k,1)&RevTempEdges(:,2)==TempEdges(k,2));
        if isempty(NoOfEdge)
            MultiplicationOfTempEdges(k)=1;
        else
            MultiplicationOfTempEdges(k)=2;
            MultiplicationOfTempEdges(NoOfEdge)=2;
        end
    end
end

ContourEdges=TempEdges(MultiplicationOfTempEdges==1,:);
NoOfContourEdges=size(ContourEdges,1);

% evaluation of the regions
NoOfRegions=1;
Regions{NoOfRegions}=ContourEdges(1,1);
RefNode=ContourEdges(1,2);
ContourEdges(1,:)=[];
while size(ContourEdges,1)>0
    IndexOfNextEdge=find(ContourEdges(:,1)==RefNode);  
    if isempty(IndexOfNextEdge)
        Regions{NoOfRegions}=[Regions{NoOfRegions} RefNode];
        if size(ContourEdges,1)>0
            NoOfRegions=NoOfRegions+1;
            Regions{NoOfRegions}=ContourEdges(1,1);
            RefNode=ContourEdges(1,2);
            ContourEdges(1,:)=[];            
        end
    else
        if size(IndexOfNextEdge,1)>1      
            IndexOfNextEdge;
            PrevNode= Regions{NoOfRegions}(end);
            TempNodes=ContourEdges(IndexOfNextEdge,2);            
            Index= FindNextNode(NodesCoord,PrevNode,RefNode,TempNodes);           
            IndexOfNextEdge=IndexOfNextEdge(Index);
	   end
        NextEdge=ContourEdges(IndexOfNextEdge,:);
        Regions{NoOfRegions}=[Regions{NoOfRegions} ContourEdges(IndexOfNextEdge,1)];
        RefNode=ContourEdges(IndexOfNextEdge,2);
        ContourEdges(IndexOfNextEdge,:)=[];
    end
end

Regions{NoOfRegions}=[Regions{NoOfRegions} RefNode];

for k = 1:NoOfRegions
    QuadrantSequence=Quadrants(Regions{k});
    dQ=QuadrantSequence(2:end)-QuadrantSequence(1:end-1);
    dQ(dQ==3)=-1;
    dQ(dQ==-3)=1;
    dQ(abs(dQ)==2)=NaN;
    q(k,1)=sum(dQ)/4;
    z(k,1)=mean(NodesCoord(Regions{k},1)+1i*NodesCoord(Regions{k},2));
end

warning('on','all');
z_root=z(q>0);
z_roots_multiplicity=q(q>0);
nrts=size(z_root,1);
[~,id] = sort(imag(z_root));     % sort roots on imaginary part
rts = z_root(id);
rtsMul = z_roots_multiplicity(id);










