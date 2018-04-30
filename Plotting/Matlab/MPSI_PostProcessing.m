close all; clear all;
%% PLOT ION DENSITY
width=750;
height=400;
set(gcf,'units','points','position',[10,10,width,height])

ax1=axes; hold on;
ax2=axes;

%set([ax1,ax2],'Position',[.17 .11 .685 .815]);
set([ax1,ax2],'Position',[.19 .11 .685 .815]);
datatmp = importdata('../../PlasmaData/Magnum-PSI/Magnum-PSI_Experiment/Magnum-PSI_Prelim_B1.41_L1.0_IonDensity.txt');
xvec=datatmp.data(1,:);
yvec=cellfun(@(x) str2double(x),datatmp.textdata(2:end));
mat = datatmp.data(2:end,:);
imagesc(ax1,xvec,yvec,mat);
imagesc(ax1,xvec,-yvec,mat);
cgray = colormap(ax1,gray(1000));
cgray = flipud(cgray);
cmap1=colormap(ax1,cgray);

%% READ PRELIMINARY DATA FILES
%ForceData = readtable('Data/Magnum-PSI_Prelim_B1.41_L1.0/3.25um/MPSI_1_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Prelim_B1.41_L1.0/3.25um/MPSI_1_hm_0.txt');

%% READ DRAFT PLOT DATA FILES
ForceData = readtable('../../Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9/3.25um/MPSI_3_fm_0.txt');
HeatData = importdata('../../Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9/3.25um/MPSI_3_hm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9/8.0um/MPSI_3_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9/8.0um/MPSI_3_hm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9/3.25um/MPSI_3_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9/3.25um/MPSI_3_hm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9/8.0um/MPSI_5_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9/8.0um/MPSI_5_hm_0.txt');

%% READ NEW DATA FILES
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/3.25um/MPSI_5_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/3.25um/MPSI_5_hm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/4.5um/MPSI_5_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/4.5um/MPSI_5_hm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/8.0um/MPSI_5_fm_0.txt');
%HeatData = importdata('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/8.0um/MPSI_5_hm_0.txt');

%% READ THE DATA

SizeOfFData=size(ForceData);
PositionData=zeros(SizeOfFData(1),3);
ForceTime=table2array(ForceData(:,1));
PositionData(:,1)=table2array(ForceData(:,2));
PositionData(:,2)=table2array(ForceData(:,3));
PositionData(:,3)=table2array(ForceData(:,4));
VelocityData=table2array(ForceData(:,7));
Temp=HeatData.data(:,2);
HeatTime=HeatData.data(:,1);
myColorMap=jet(256);
n=length(Temp(:));
cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
%Tempmax=max(Temp(:));

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
fprintf('Tempmax = %d',Tempmax);

%% PLOT LINES
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

for k = 1 : (size(PositionData(:,1))-6)
	%plot( ax2, Z([k k+1],1),X([k k+1],1),'Color', theLineColor(k,:), 'LineWidth',8.0);
	line( ax2,[Z(k,1) Z(k+5,2)], [X(k,1) X(k+5,2)],'Color', theLineColor(k,:), 'LineWidth',2.5);
	if ~(mod(k,1000))
		hold on;
		plotm1 = plot(ax2, Z(k,1),X(k,1),'k*','markers',12,'LineWidth',3.0);
		%plotm2 = plot(ax2, Z(k,1),X(k,1),'ko','markers',12,'LineWidth',3.0);
		%plotm3 = plot(ax2, Z(k,1),X(k,1),'kx','markers',12,'LineWidth',3.0);
		%plotm4 = plot(ax2, Z(k,1),X(k,1),'ks','markers',12,'LineWidth',3.0);
		%plotm5 = plot(ax2, Z(k,1),X(k,1),'kd','markers',12,'LineWidth',3.0);
		%plotm6 = plot(ax2, Z(k,1),X(k,1),'k+','markers',12,'LineWidth',3.0);
	end
end
%% design
%legend([plotm],{'B0.1 R8.0\mum'},'location','NorthWest');
% link axis together
linkaxes([ax1,ax2]);
ax2.Visible='off';
ax2.XTick=[];
ax2.YTick=[];


% Sort out color bars

cmap2=colormap(ax2,'jet');

xt = get(ax1, 'XTick');
yt = get(ax1, 'YTick');
set(ax1, 'FontSize', 20);

xlabel('Longitudinal Position (m)','FontSize',32);
%ylabel('Vertical Position (m)','FontSize',32);
grid on;
%cb1=colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb1=colorbar(ax1,'Position',[.06 .11 .0675 .815]);
cb1.Label.String = 'Ion Density (m^-^3)';
cb1.FontSize = 20;
%cb2=colorbar(ax2,'Position',[.88 .11 .0675 .815]);
cb2=colorbar(ax2,'Position',[.82 .11 .0675 .815]);
cb2.Ticks = linspace(0, 1, 7) ; %Create 8 ticks from zero to 1
cb2.TickLabels = num2cell(round(min(Temp)):450:Tempmax);    %Replace the labels of these 8 ticks with the numbers 1 to 8
cb2.Label.String = 'Temperature (K)';
cb2.FontSize = 20;
xlim([0.1 0.65])
