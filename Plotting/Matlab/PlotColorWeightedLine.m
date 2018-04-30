function [  ] = PlotColorWeightedLine( ForceData, HeatData, Tempmax, ax2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

SizeOfFData=size(ForceData);
PositionData=zeros(SizeOfFData(1),3);
ForceTime=table2array(ForceData(:,1));
PositionData(:,1)=table2array(ForceData(:,2));
PositionData(:,2)=table2array(ForceData(:,3));
PositionData(:,3)=table2array(ForceData(:,4));
Temp=HeatData.data(:,2);
HeatTime=HeatData.data(:,1);
myColorMap=jet(256);
n=length(Temp(:));
cd = [uint8(jet(n)*255) uint8(ones(n,1))].';


%[X,Y,Z]=cylinder(0.15);
%[X2,Y2,Z2]=cylinder(0.25);
%s1=surf(X,Y,Z);
%s2=surf(X2,Y2,Z2);
RadialPos=zeros(length(PositionData(:,1)),2);
Theta=zeros(length(PositionData(:,1)),2);
X=zeros(length(PositionData(:,1)),2);
Y=zeros(length(PositionData(:,1)),2);
Z=zeros(length(PositionData(:,1)),2);
theIntensity=zeros(length(PositionData(:,1)));
theLineColor=zeros(length(PositionData(:,1)),3);
j=1;

for k = 1 : size(PositionData(:,1))-1
	
	RadialPos(k,1) = PositionData(k,1);
	RadialPos(k,2) = PositionData(k+1,1);
	Theta(k,1) = PositionData(k,2);
	Theta(k,2) = PositionData(k+1,2);

	X(k,1) = PositionData(k,1)*cos(PositionData(k,2));
	X(k,2) = PositionData(k+1,1)*cos(PositionData(k+1,2));
	Y(k,1) = PositionData(k,1)*sin(PositionData(k,2));	
	Y(k,2) = PositionData(k+1,1)*sin(PositionData(k+1,2));
	Z(k,1) = PositionData(k,3);
	Z(k,2) = PositionData(k+1,3);
	if( ForceTime(k) >= HeatTime(j) )
		j = j +1;
	end
	
	theIntensity(k)=Temp(j)/Tempmax;
	theLineColor(k,:)=myColorMap(round(theIntensity(k)*256),:);
end

for k = 1 : size(PositionData(:,1))-1
	%plot( ax2, Z([k k+1],1),X([k k+1],1),'Color', theLineColor(k,:), 'LineWidth',8.0);
	line( ax2,[Z(k,1) Z(k,2)], [X(k,1) X(k,2)],'Color', theLineColor(k,:));%, 'LineWidth',1.0);
end

end

