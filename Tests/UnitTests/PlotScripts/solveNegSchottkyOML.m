function [ minimized_x ] = solveNegSchottkyOML( Ti, Te, Td, ni, ne, Dsec, init_guess )
%SOLVEPOSSCHOTTKYOML Summary of this function goes here
%   Detailed explanation goes here
        echarge = 1.6e-19;
        Me = 9.1e-31;
        Mp = 1.66e-27;
        Kb = 1.38e-23;
		Richardson = 1.2e6;
        
		WorkFunction = 3.4*echarge;     
		Ve = sqrt(Kb*Te/Me);
		Vi = sqrt(Kb*Ti/Mp);

		a = ne*Ve*(Dsec-1);
		b = ni*Vi*(Te/Ti); 
		c = ni*Vi;
		d = 4*Richardson*(Td^2)/echarge;
		f = WorkFunction/(Kb*Td);

		%fx      = @(guess) (a.*exp(-guess)+b.*guess+C+d.*exp(guess-f)).^2;
		%fx      = @(guess) (a.*exp(-guess)+b.*guess+c+d.*exp(guess-f)).^2;
		fx      = @(guess) ((a./c).*exp(-guess)+(b./c).*guess+1+(d./c).*exp(guess-f)).^2;
		% debug code:
%		fprintf('\na = %1.2e\n',a);
%		fprintf('\nb = %1.2e\n',b);
%		fprintf('\nc = %1.2e\n',c);
%		fprintf('\nd = %1.2e\n',d);
%		fprintf('\nf = %1.2e\n',f);
		% init_guess = -0.1;
 		j=1;
% 		for init_guess = 0.1:0.1:30
%         opt = optimset('TolX',1e-15,'MaxFunEvals',1e6);
        [minimized_x(j),fval(j)] = fminsearch(fx,init_guess);
% 		j=j+1;
% 		end
% 		min(minimized_x)
% 		min(fval)
		
end