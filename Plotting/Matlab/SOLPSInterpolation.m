clear all;

datatmp = importdata('../../PlasmaData/JET/Data_2.0e19_2.2MW_noheader.txt');
R=datatmp(:,3);
Z=datatmp(:,4);
ToroidalB=datatmp(:,5);
PoloidalB=datatmp(:,6);
CellVolume=datatmp(:,7);
IonTemperature=datatmp(:,8);
ElectronTemperature=datatmp(:,9);
MolecularDensity=datatmp(:,10);
AtomicDensity=datatmp(:,11);
IonDensity=datatmp(:,12);
ElectronDensity=datatmp(:,13);
ParallelIonVelocity=datatmp(:,14);
ParallelElectronVelocity=datatmp(:,15);
PerpendicularIonVelocity=datatmp(:,16);
PerpendicularElectronVelocity=datatmp(:,17);

Xq=[1.5:0.01:4.0];
Yq=[-2.0:0.01:2.0];

ToroidalB_new=zeros(length(Xq),length(Yq));
PoloidalB_new=zeros(length(Xq),length(Yq));
CellVolume_new=zeros(length(Xq),length(Yq));
IonTemperature_new=zeros(length(Xq),length(Yq));
ElectronTemperature_new=zeros(length(Xq),length(Yq));
MolecularDensity_new=zeros(length(Xq),length(Yq));
AtomicDensity_new=zeros(length(Xq),length(Yq));
IonDensity_new=zeros(length(Xq),length(Yq));
ElectronDensity_new=zeros(length(Xq),length(Yq));
ParallelIonVelocity_new=zeros(length(Xq),length(Yq));
ParallelElectronVelocity_new=zeros(length(Xq),length(Yq));
PerpendicularIonVelocity_new=zeros(length(Xq),length(Yq));
PerpendicularElectronVelocity_new=zeros(length(Xq),length(Yq));


Tolerance=0.1;

Weightings = zeros(length(Xq),length(Yq),length(R));
indices = zeros(length(Xq),length(Yq),length(R));

for i = 1:length(Xq)
	for j = 1:length(Yq)
		sprintf('X : %i, Y : %i',i,j)
		num_in_range = 0;
		for k = 1:length(R)
			%IrregularGridPosition=sqrt(R(k)*R(k)+Z(k)*Z(k));
			%RegularGridPosition=sqrt(Xq(i)*Xq(i)+Yq(i)*Yq(i);
			rdisp = (Xq(i)-R(k));
			zdisp = (Yq(j)-Z(k));
			Seperation=sqrt(rdisp*rdisp+zdisp*zdisp);
			if ( Seperation < Tolerance )
				num_in_range = num_in_range + 1;
				Weightings(i,j,num_in_range)=Seperation;
				indices(i,j,num_in_range)=k;
			end
		end
		for p = 1:length(Weightings(i,j,:))
			if( length(Weightings(i,j,:)) > 1 )
				Weight = Weightings(i,j,p);
				ToroidalB_new(i,j)=ToroidalB_new(i,j)+Weight*ToroidalB(k)/p;
				PoloidalB_new(i,j)=PoloidalB_new(i,j)+Weight*PoloidalB(k)/p;
				CellVolume_new(i,j)=CellVolume_new(i,j)+Weight*CellVolume(k)/p;
				IonTemperature_new(i,j)=IonTemperature_new(i,j)+Weight*IonTemperature(k)/p;
				ElectronTemperature_new(i,j)=ElectronTemperature_new(i,j)+Weight*ElectronTemperature(k)/p;
				MolecularDensity_new(i,j)=MolecularDensity_new(i,j)+Weight*MolecularDensity(k)/p;
				AtomicDensity_new(i,j)=AtomicDensity_new(i,j)+Weight*AtomicDensity(k)/p;
				IonDensity_new(i,j)=IonDensity_new(i,j)+Weight*IonDensity(k)/p;
				ElectronDensity_new(i,j)=ElectronDensity_new(i,j)+Weight*ElectronDensity(k)/p;
				ParallelIonVelocity_new(i,j)=ParallelIonVelocity_new(i,j)+Weight*ParallelIonVelocity(k)/p;
				ParallelElectronVelocity_new(i,j)=ParallelElectronVelocity_new(i,j)+Weight*ParallelElectronVelocity(k)/p;
				PerpendicularIonVelocity_new(i,j)=PerpendicularIonVelocity_new(i,j)+Weight*PerpendicularIonVelocity(k)/p;
				PerpendicularElectronVelocity_new(i,j)=PerpendicularElectronVelocity_new(i,j)+Weight*PerpendicularElectronVelocity(k)/p;
			end
		end
	end
end
 

% [Sorted_R,R_indices]=sort(R,'ascend');
% [Sorted_Z,Z_indices]=sort(Z,'ascend');
% RPrime = unique(R);
% ZPrime = unique(Z);
% 
% [RMesh,ZMesh] = meshgrid(RPrime,ZPrime');
% ToroidalBMatrix = NaN(size(RMesh)); 
% PoloidalBMatrix = NaN(size(RMesh)); 
% for i = 1:length(ToroidalB)
% 	idR = find( RMesh(1,:) == R(i) );
% 	idZ = find( ZMesh(:,1) == Z(i) );
%  	
%  	ToroidalBMatrix(idR,idZ) = ToroidalB(i);
%  	
%  	PoloidalBMatrix(idR,idZ) = PoloidalB(i);
% end
% ToroidalBMatrix_Painted = inpaint_nans(ToroidalBMatrix);
% 
% ToroidalBMatrix_Painted = inpaint_nans(PoloidalBMatrix);
% figure(1);
% imagesc(ToroidalBMatrix)
% figure(2)
% imagesc(ToroidalBMatrix_Painted)
% 
% ToroidalB_new=interpn(Sorted_R,Sorted_Z,ToroidalB,Xq,Yq);











%PoloidalB_new=interp2(R,Z,PoloidalB,Xq,Yq);