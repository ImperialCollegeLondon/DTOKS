k(1) = 2.5;

% Two different initial guesses for the minimisation search algorithm
m(1) = -0.1;
p(1) = -0.001;

for i = 300 : 5555
	k(i-298) = solveNegSchottkyOML( 10e5, 10e5, i, 1e18, 1e18, 0.01, k(i-299) );
	m(i-298) = solvePosSchottkyOML( 10e5, 10e5, i, 1e18, 1e18, 0.01, m(i-299) );
	p(i-298) = solvePosSchottkyOML( 10e5, 10e5, i, 1e18, 1e18, 0.01, p(i-299) );
end


plot(k(k>0));
hold on
% This arbitary offset was implemented to make the two solutions match.
% It's not clear why it is necessary or how it enters into the problem.
m=m+1.86;
p=p+1.86;
plot(m);

% For some strange reason, the average of the two offset branches is the
% correct solution that maps onto the solution of the negative equation.
% I can't explain this further and am leaving this problem here without a
% full solution.
plot((p+m)/2);