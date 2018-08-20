function [ minimized_x ] = solvePosSchottkyOML( Ti, Te, Td, ni, ne, Dsec, init_guess )
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

        Esee = 3*echarge;
        K = Kb*Te/Esee;
        A = ne*Ve*Dsec;
        B = (Te*Dsec)/Ti;
        C = -ne*Ve;
        D = -Te/Ti;
        H = ni*Vi; 
        I = 4*Richardson*(Td^2)/echarge;
        f = WorkFunction/(Kb*Td);

        %fx      = @(guess) C+A*exp(K*guess)+B*guess*exp(K*guess)+D*guess+H*exp(-guess)+I*exp(guess-f);
		fx      = @(guess) (1+(A./C)*exp(K*guess)+(B./C)*guess*exp(K*guess)+(D./C)*guess+(H./C)*exp(-guess)+(I./C)*exp(guess-f)).^2;
%		fprintf('\nA = %1.2e\n',A);
%		fprintf('\nB = %1.2e\n',B);
%		fprintf('\nC = %1.2e\n',C);
%		fprintf('\nD = %1.2e\n',D);
%		fprintf('\nH = %1.2e\n',H);
%		fprintf('\nI = %1.2e\n',I);
        % init_guess = -0.1;
        % opt = ('TolX',1e-15,'MaxFunEvals',1e6);
        
        minimized_x = fminsearch(fx,init_guess);
        
end

