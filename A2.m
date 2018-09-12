%A2
clear
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
gamma = 1;
Smax=4*K;

t_num=2000; %num of time-step
M=100;

t=linspace(0,T,t_num);
s=linspace(0,Smax,M);
delta_t=T/t_num;
delta_s=Smax/M;

last_v=max(s-K,0);%?????
%% Euler backwards
for n=t_num:-1:t_num-1
    v(1)=0;
    v(M)=Smax-K*exp(-r*(T-t(n-1)));
    for j=2:M-1
        v(j)=last_v(j)+r*s(j)*delta_t/(2*delta_s)*(last_v(j+1)-last_v(j-1))...
             +sigma^2*0.5*(s(j))^2*(delta_t/(delta_s^2))*(last_v(j+1)-...
             2*last_v(j)+last_v(j-1))-delta_t*r*last_v(j);
    end
    last_v=v;
end
%% Crank-Nicholson 



%%  compare with exact solution
for i=1:M
    exact(i) = bsexact(sigma, r, K, 0, s(i)); %?????? 
end
%% plots
plot(s,last_v,'r*')
hold on
plot(s,exact,'b')

figure
plot(s,last_v-exact,'r')
