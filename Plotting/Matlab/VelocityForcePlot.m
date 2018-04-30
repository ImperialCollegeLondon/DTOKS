ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9/3.25um/MPSI_3_fm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9/8.0um/MPSI_3_fm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9/3.25um/MPSI_3_fm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9/8.0um/MPSI_5_fm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/3.25um/MPSI_5_fm_0.txt');
%ForceData = readtable('Data/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9/8.0um/MPSI_5_fm_0.txt');

SizeOfFData=size(ForceData);
PositionData=zeros(SizeOfFData(1),3);

PositionData(:,1)=table2array(ForceData(:,2));
PositionData(:,2)=table2array(ForceData(:,3));
PositionData(:,3)=table2array(ForceData(:,4));
VelocityDatax=table2array(ForceData(:,5));
VelocityDatay=table2array(ForceData(:,6));
VelocityDataz=table2array(ForceData(:,7));
Lorentzz=table2array(ForceData(:,17));
IonDragz=table2array(ForceData(:,20));
ForceTime=table2array(ForceData(:,1));

X = PositionData(:,1).*cos(PositionData(:,2));
figure(1);
hold on;
%plotm11 = plot(ForceTime,IonDragz/50,'k','LineWidth',3.0);
%plotm12 = plot(ForceTime,IonDragz/50,'r','LineWidth',3.0);
%plotm13 = plot(ForceTime,IonDragz/50,'b','LineWidth',3.0);
%plotm14 = plot(ForceTime,IonDragz/50,'g','LineWidth',3.0);
for k = 1 : (size(PositionData(:,1))-6)
	
	if ~(mod(k,1000))
		hold on;
		
		plotm1 = plot(ForceTime(k),VelocityDataz(k),'k*','markers',12,'LineWidth',3.0);
		%plotm2 = plot(ForceTime(k),VelocityDataz(k),'ko','markers',12,'LineWidth',3.0);
		%plotm3 = plot(ForceTime(k),VelocityDataz(k),'kx','markers',12,'LineWidth',3.0);
		%plotm4 = plot(ForceTime(k),VelocityDataz(k),'ks','markers',12,'LineWidth',3.0);
		%plotm5 = plot(ForceTime(k),VelocityDataz(k),'kd','markers',12,'LineWidth',3.0);
		%plotm6 = plot(ForceTime(k),VelocityDataz(k),'k+','markers',12,'LineWidth',3.0);
		%plotm41 = plot(ForceTime(k),IonDragz(k)/50,'gs','markers',12,'LineWidth',3.0);
		%marker_order=marker_order+1;
	end
end
%legend([plotm1 plotm2 plotm3 plotm4 plotm5 plotm6],{'B=0.1 R=3.25\mum','B=0.1 R=8.0\mum','B=0.2 R=3.25\mum','B=0.2 R=8.0\mum','B=0.4 R=3.25\mum','B=0.4 R=8.0\mum'},'location','NorthWest');
xlabel('Time (s)')
ylabel('Horizontal Velocity (ms^-^1)')
grid on