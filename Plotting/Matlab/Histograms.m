close all;
clear all;

%% Data Input
NumFiles=9;
NumberOfLines=1000;
fid=zeros(NumFiles,1);
FinalPosz=zeros(1000,NumFiles);
filenames={
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9_3.25um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9_4.5um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9_8.0um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9_3.25um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9_4.5um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9_8.0um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9_3.25um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9_4.5um.txt', ...
	'Data/FinalPos_Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9_8.0um.txt'};


for i=1:NumFiles
	fid(i)=fopen(filenames{i});
	TextScanData=textscan(fid(i),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter',' \t\n');
	FinalPosz(:,i)=cell2mat(TextScanData(4))-0.15;
end
fclose('all');




%% Error Analysis
a1=sort(FinalPosz(:,1));
a2=sort(FinalPosz(:,2));
a3=sort(FinalPosz(:,3));
a4=sort(FinalPosz(:,4));
a5=sort(FinalPosz(:,5));
a6=sort(FinalPosz(:,6));
a7=sort(FinalPosz(:,7));
a8=sort(FinalPosz(:,8));
a9=sort(FinalPosz(:,9));

[f0, x0] = hist(FinalPosz(:,1)); % Create histogram from a normal distribution.
[f1, x1] = hist(FinalPosz(:,2)); % Create histogram from a normal distribution.
[f2, x2] = hist(FinalPosz(:,3)); % Create histogram from a normal distribution.
[f3, x3] = hist(FinalPosz(:,4)); % Create histogram from a normal distribution.
[f4, x4] = hist(FinalPosz(:,5)); % Create histogram from a normal distribution.
[f5, x5] = hist(FinalPosz(:,6)); % Create histogram from a normal distribution.
[f6, x6] = hist(FinalPosz(:,7)); % Create histogram from a normal distribution.
[f7, x7] = hist(FinalPosz(:,8)); % Create histogram from a normal distribution.
[f8, x8] = hist(FinalPosz(:,9)); % Create histogram from a normal distribution.

err1=zeros(10,1);
err2=err1;
err3=err1;
err4=err1;
err5=err1;
err6=err1;
err7=err1;
err8=err1;
err9=err1;
for i = 1:10
	if( i == 1 )
		err1(i)=sqrt(sum((a1(1:f0(1))-x0(i)).^2)/f0(1));
		err2(i)=sqrt(sum((a2(1:f1(1))-x1(i)).^2)/f1(1));
		err3(i)=sqrt(sum((a3(1:f2(1))-x2(i)).^2)/f2(1));
		err4(i)=sqrt(sum((a4(1:f3(1))-x3(i)).^2)/f3(1));
		err5(i)=sqrt(sum((a5(1:f4(1))-x4(i)).^2)/f4(1));
		err6(i)=sqrt(sum((a6(1:f5(1))-x5(i)).^2)/f5(1));
		err7(i)=sqrt(sum((a7(1:f6(1))-x6(i)).^2)/f6(1));
		err8(i)=sqrt(sum((a8(1:f7(1))-x7(i)).^2)/f7(1));
		err9(i)=sqrt(sum((a9(1:f8(1))-x8(i)).^2)/f8(1));
	else
		err1(i)=sqrt(sum((a1(sum(f0(1:i-1)):sum(f0(1:i)))-x0(i)).^2)/f0(i));
		err2(i)=sqrt(sum((a2(sum(f1(1:i-1)):sum(f1(1:i)))-x1(i)).^2)/f1(i));
		err3(i)=sqrt(sum((a3(sum(f2(1:i-1)):sum(f2(1:i)))-x2(i)).^2)/f2(i));
		err4(i)=sqrt(sum((a4(sum(f3(1:i-1)):sum(f3(1:i)))-x3(i)).^2)/f3(i));
		err5(i)=sqrt(sum((a5(sum(f4(1:i-1)):sum(f4(1:i)))-x4(i)).^2)/f4(i));
		err6(i)=sqrt(sum((a6(sum(f5(1:i-1)):sum(f5(1:i)))-x5(i)).^2)/f5(i));
		err7(i)=sqrt(sum((a7(sum(f6(1:i-1)):sum(f6(1:i)))-x6(i)).^2)/f6(i));
		err8(i)=sqrt(sum((a8(sum(f7(1:i-1)):sum(f7(1:i)))-x7(i)).^2)/f7(i));
		err9(i)=sqrt(sum((a9(sum(f8(1:i-1)):sum(f8(1:i)))-x8(i)).^2)/f8(i));
	end
end

%% Figure 1: Small Magnetic Field, many radii
% METHOD 1: DIVIDE BY SUM

%g = 1 / sqrt(2 * pi) * exp(-0.5 * x .^ 2); % pdf of the normal distribution
figure(1);
p1=bar(x0, f0 / sum(f0)); hold on
p2=bar(x1, f1 / sum(f1)); hold on
p3=bar(x2, f2 / sum(f2)); hold on
set(p1,'FaceColor','red','EdgeColor','black','FaceAlpha',.8);
set(p2,'FaceColor','blue','EdgeColor','black','FaceAlpha',.8);
set(p3,'FaceColor','green','EdgeColor','black','FaceAlpha',.8);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 30)
title('Final Positions B=0.1T','FontSize',36);
xlabel('Longitudinal Position (m)','FontSize',32);
ylabel('Normalised Frequency (arb)','FontSize',32);
grid on;
legend('16\mum','9\mum','6.5\mum');


%% Figure 2: Large Magnetic Field, many radii
figure(2)
p4=bar(x3, f3 / sum(f3)); hold on
%p5=bar(x4, f4 / sum(f4)); hold on
p6=bar(x5, f5 / sum(f5)); hold on

%eb1=errorbar(x0,f0/sum(f0),err1,'.k'); set(eb1, 'MarkerSize', 10);
%eb2=errorbar(x1,f1/sum(f1),err2,'.k'); set(eb2, 'MarkerSize', 10);
%eb3=errorbar(x2,f2/sum(f2),err3,'.k'); set(eb3, 'MarkerSize', 10);
eb4=errorbar(x3,f3/sum(f3),err4,'.k'); set(eb4, 'MarkerSize', 10);
%eb5=errorbar(x4,f4/sum(f4),err5,'.k'); set(eb5, 'MarkerSize', 10);
eb6=errorbar(x5,f5/sum(f5),err6,'.k'); set(eb6, 'MarkerSize', 10);

set(p4,'FaceColor','red','EdgeColor','black','FaceAlpha',.8);
%set(p5,'FaceColor','green','EdgeColor','black','FaceAlpha',.8);
set(p6,'FaceColor','blue','EdgeColor','black','FaceAlpha',.8);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 30)
title('Final Positions B=0.2T','FontSize',36);
xlabel('Longitudinal Position (m)','FontSize',32);
ylabel('Normalised Frequency (arb)','FontSize',32);
grid on;
legend('6.5\mum','16\mum');


%% Figure 3: Larger Magnetic Field, many radii
figure(3)
%p7=bar(x6, f6 / sum(f6)); hold on
%p8=bar(x7, f7 / sum(f7)); hold on
%p9=bar(x8, f8 / sum(f8)); hold on
%set(p7,'FaceColor','red','EdgeColor','black','FaceAlpha',.8,'BarWidth',0.5);
%set(p8,'FaceColor','blue','EdgeColor','black','FaceAlpha',.8,'BarWidth',0.5);
%set(p9,'FaceColor','green','EdgeColor','black','FaceAlpha',.8,'BarWidth',0.8);
p7=plot(x6, f6 / sum(f6),'rx','markers',12,'LineWidth',3.0); hold on;
p8=plot(x7, f7 / sum(f7),'bx','markers',12,'LineWidth',3.0); hold on;
p9=plot(x8, f8 / sum(f8),'gx','markers',12,'LineWidth',3.0); hold on;
eb7=errorbar(x6,f6/sum(f6),err7,'rx'); set(eb7,'markers',12, 'MarkerSize', 10);
eb8=errorbar(x7,f7/sum(f7),err8,'bx'); set(eb8,'markers',12, 'MarkerSize', 10);
eb9=errorbar(x8,f8/sum(f8),err9,'gx'); set(eb9,'markers',12, 'MarkerSize', 10);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 30)
title('Final Positions B=0.4T','FontSize',36);
xlabel('Longitudinal Position (m)','FontSize',32);
ylabel('Normalised Frequency (arb)','FontSize',32);
grid on;
legend('3.25\mum','4.5\mum','8.0\mum');

%% Figure 4: Small radii, Many magnetic fields
figure(4)
%p1=bar(x0, f0 / sum(f0)); hold on
%p4=bar(x3, f3 / sum(f3)); hold on
%p7=bar(x6, f6 / sum(f6)); hold on
%p1=plot(x0, f0 / sum(f0),'rx','markers',12,'LineWidth',3.0); hold on;
p4=plot(x3, f3 / sum(f3),'bx','markers',12,'LineWidth',3.0); hold on;
p7=plot(x6, f6 / sum(f6),'gx','markers',12,'LineWidth',3.0); hold on;
%eb1=errorbar(x0,f0/sum(f0),err1,'.k'); set(eb1, 'MarkerSize', 10);
eb4=errorbar(x3,f3/sum(f3),err4,'.k'); set(eb4, 'MarkerSize', 10);
eb7=errorbar(x6,f6/sum(f6),err7,'.k'); set(eb7, 'MarkerSize', 10);

%set(p1,'FaceColor','red','EdgeColor','black','FaceAlpha',.8,'BarWidth',1.0);
%set(p4,'FaceColor','blue','EdgeColor','black','FaceAlpha',.8,'BarWidth',1.0);
%set(p7,'FaceColor','green','EdgeColor','black','FaceAlpha',.8,'BarWidth',0.5);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 30)
title('Final Positions r=3.25\mum','FontSize',36);
xlabel('Longitudinal Position (m)','FontSize',32);
ylabel('Normalised Frequency (arb)','FontSize',32);
grid on;
legend('0.2T','0.4T'); %'0.1T',