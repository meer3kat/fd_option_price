%A2
clear
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
gamma = 1;
Smax = 4 * K;

t_num = 100; %num of time-step
M = 61;

t = linspace(0,T,t_num);
s = linspace(0,Smax,M);
delta_t = t(2)-t(1);
delta_s = s(2)-s(1);

last_v = max(s-K,0); 
%% Euler explicit
for n=t_num:-1:2
    v(1) = 0;
    v(M) = Smax-K*exp(-r*(T-t(n-1)));
    for j=2:M-1
        v(j)=last_v(j)+r*s(j)*delta_t/(2*delta_s)*(last_v(j+1)-last_v(j-1))...
             +sigma^2*0.5*(s(j))^2*(delta_t/(delta_s^2))*(last_v(j+1)-...
             2*last_v(j)+last_v(j-1))-delta_t*r*last_v(j);
    end
    last_v=v;
end
%% Euler implicit


v_i = max(s-K,0); 

%v_e = max(s-K,0);%?



for n=t_num:-1:2
    
    vi(1) = 0;
    %ve(1,M)=Smax-K*exp(-r*(T-t(n-1)));
    vi(M)=Smax-K*exp(-r*(T-t(n-1)));
    
    
    ci = r*s/(2*delta_s) + 0.5*sigma^2*s.^2/(delta_s*delta_s);
    bi = -sigma^2*s.^2 /(delta_s*delta_s) - r-1/delta_t;
    ai = -r*s/(2*delta_s) + 0.5*sigma^2*s.^2/(delta_s*delta_s);
    
    D = diag(ai(3:end-1),-1)+ diag(bi(2:end-1)) + diag(ci(2:end-2),1); %
    %bm = v_i(2:end-1);
    %bm = bm(:); 
    %bm(1) = bm(1) - ai(1)*vi(1);
    %bm(end) = bm(end) - ci(end)*vi(M);
    bm=-1/delta_t*v_i(2:end-1);
    bm(1) = bm(1) - ai(2)*vi(1);
    bm(end) = bm(end) - ci(end-1)*vi(M);
    vi(2:end-1) = D\bm';
    v_i = vi;
    
    
    %D = diag(ai(2:end-1),-1) + diag(bi(1:end-1)) + diag(ci(1:end-2),1);
    %B = eye(size(D)) + delta_t*D;
    %B = B(2:end-1,2:end-1);
    %A = eye(size(D)) - delta_t*D;
    %A = A(2:end-1,2:end-1);

end
 


%%  compare with exact solution
for i=1:M
    exact(i) = bsexact(sigma, r, K, T, s(i)); 
end
%% plots
%plot(s,last_v,'r*')
plot(s,vi,'gd')
hold on
plot(s,exact,'b')


 figure
 plot(s,abs(vi-exact),'r')
 hold on
 plot(s,abs(v-exact),'b')