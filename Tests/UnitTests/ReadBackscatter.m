function [ Time, eTemp, iTemp, Potential, RE, RN ] = ReadBackscatter( filename )
% Function to read in data produced by HeatingMatter program depending on
% which heating models are activated and which aren't

Import = importdata(filename);

Time = Import(:,1);
eTemp = Import(:,2);
iTemp = Import(:,3);
Potential = Import(:,4);
RE = Import(:,5);
RN = Import(:,6);

end