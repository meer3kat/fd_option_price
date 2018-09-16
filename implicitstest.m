clear all;
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
g = 0.9;
Smax = 4 * K;
tn = 101;
sn = 61;

[s,vj] = fdimplicit(K,r,sigma,T,g,tn,sn);

%  compare with exact solution
for i=1:length(s)
    exact(i) = bsexact(sigma, r, K, T, s(i)); %?????? 
end
%% plots
plot(s,vj,'r*')
hold on
plot(s,exact,'b')
%plot(s,v_tn,'gd')

err = abs(exact - vj');
figure
plot(s,err,'r')
