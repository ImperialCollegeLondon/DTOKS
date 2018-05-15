clear all;
close all;
grid on
xlabel('Major Radius (m)')
ylabel('Toroidal Angle (m)');
%ylabel('Z (m)');
zlabel('Height (m)');
set(get(gca,'zlabel'),'rotation',0)
%view([0.3 -0.6 0.1]);
view([-0.5 -0.5 0.2]);
%dirname='Data/BreakupDataForPaper/new/ITER_10um_Bdir_W_4/';
dirname='/home/ls5115/Software/DTOKS-U/Data/';
filenum=33;
ForceData = readtable(strcat(dirname,'MAST_fm_0.txt'));
HeatData = readtable(strcat(dirname,'MAST_hm_0.txt'));

Density= 19600;	% Density of Tungsten in kg/m^3

SizeOfFData=size(ForceData);
PositionData=zeros(SizeOfFData(1),3);

M0=table2array(HeatData(1,3));
R0=(M0*3.0/(Density*4.0*pi))^(1.0/3.0);

TimeData=table2array(ForceData(:,1));
Times=zeros(1,1520);
Times(1)=TimeData(end)-TimeData(1);
for i = 1:filenum
	sprintf('%i',i)
	ForceData = readtable(strcat(dirname,sprintf('breakup_fm_%i.txt',i)));
	TimeData=table2array(ForceData(:,1));
	Times(i+1)=TimeData(end)-TimeData(1);
end

%%
ForceData = readtable(strcat(dirname,'MAST_fm_0.txt'));
HeatData = readtable(strcat(dirname,'MAST_hm_0.txt'));

PositionData(:,1)=table2array(ForceData(:,2));
PositionData(:,2)=table2array(ForceData(:,3));
PositionData(:,3)=table2array(ForceData(:,4));
Mass=table2array(HeatData(1,3));
Radius=(Mass*3.0/(Density*4.0*pi))^(1.0/3.0);
myColorMap=jet(256);
theIntensity=Radius/R0;
theLineColor=myColorMap(ceil(theIntensity*256),:);
figure(1);
plot3(PositionData(:,1),PositionData(:,2),PositionData(:,3),'Color', theLineColor,'LineWidth',4);
%plot3(cos(PositionData(:,2)).*PositionData(:,1),sin(PositionData(:,2)).*PositionData(:,1),PositionData(:,3),'Color', theLineColor,'LineWidth',4);
hold on;


%figure(2);
%plot(PositionData(:,1),PositionData(:,3),'Color', theLineColor,'LineWidth',4);
%hold on;
for i = 1:filenum
	sprintf('%i',i)
	ForceData = readtable(strcat(dirname,sprintf('breakup_fm_%i.txt',i)));
	HeatData = readtable(strcat(dirname,sprintf('breakup_hm_%i.txt',i)));
	Mass=table2array(HeatData(1,3));
	Radius=(Mass*3.0/(Density*4.0*pi))^(1.0/3.0);
	if( Radius > 1.5e-6)
	SizeOfFData=size(ForceData);
	PositionData=zeros(SizeOfFData(1),3);
	
	PositionData(:,1)=table2array(ForceData(:,2));
	PositionData(:,2)=table2array(ForceData(:,3));
	PositionData(:,3)=table2array(ForceData(:,4));
	
	
	theIntensity=Radius/R0;
	theLineColor=myColorMap(ceil(theIntensity*256),:);
	%figure(1);
	
	plot3(PositionData(:,1),PositionData(:,2),PositionData(:,3),'Color', theLineColor,'LineWidth',4);
	%plot3(cos(PositionData(:,2)).*PositionData(:,1),sin(PositionData(:,2)).*PositionData(:,1),PositionData(:,3),'Color', theLineColor,'LineWidth',4);
	hold on;
	%figure(2);
	
	%plot(PositionData(:,1),PositionData(:,3),'Color', theLineColor,'LineWidth',4);
	%hold on;
	end
end
grid on
title('JET','fontsize',60);
xlabel('Major Radius (m)','fontsize',30);
ylabel('Toroidal Angle (rad)','fontsize',30);
%ylabel('Z (m)');
colormap jet;
hz=zlabel('Height (m)','fontsize',30);
set(get(gca,'zlabel'),'rotation',0)
view([-0.5 -0.5 0.2]);
%TickLabels_Are=[1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5,1e-4]
%cb2=colorbar('Position',[.90 .26 .02 .615]);
%cb2.Ticks = linspace(0, 1, 10) ;
%cb2.TickLabels = num2cell(TickLabels_Are)
%hold off;
