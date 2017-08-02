k(1) = 2.5;
p(1) = -0.01;
for i = 300 : 5555
	k(i-298) = solveNegSchottkyOML( 10e5, 10e5, i, 1e18, 1e18, 0.01, k(i-299) );
	p(i-298) = solvePosSchottkyOML( 10e5, 10e5, i, 1e18, 1e18, 0.01, p(i-299) );
end


plot(k(k>0));
hold on
p=p+1.86;
plot(p);