function [mu, tEnd] = GPF_UNGM(x,y)

%intial parameters
%x=0.1; %initial state
x_N=1;  %Noise variance at measurement update
x_R=1;  %Noise variance at time update
N=100;  %No. of states
M=100;   %No. of particles
%x_p_update=zeros(1,M);
%y_update=zeros(1,M);
alpha=0.5;
beta=25;gamma=8;
v=1;%variance around the initial state
x_p=zeros(1,M);
p_w = zeros(1,M);
mu = zeros(1,N);

tStart = tic;
for i=1:M
    x_p(i)=x(1)+sqrt(1)*randn;
end

for n=1:N 
    
    for j=1:M
        %updating the random variable we get
        x_p_update(j)= alpha*x_p(j) + beta*(x_p(j)/(1+x_p(j)^2)) + gamma*cos(1.2*(n-1))+sqrt(x_N)*randn;
        
        %finding the observations from the 
       % y_update(j)=x_p_update(j)^2/20+sqrt(x_R)*randn;
        
        
        p_w(j)=((1/sqrt(2*pi*x_R))*exp(-(y(n+1)-(x_p_update(j)^2/20))^2/(2*x_R)));
    end


p_w=p_w./sum(p_w);%normalizing the weights

mu_n = p_w * x_p_update';

% %Variance using the above weighted sample
var_n=0;
for k=1:M
    var_n=var_n+p_w(k)*(x_p_update(k)-mu_n)*(x_p_update(k)-mu_n);
end
% % 

for l=1:M  
    x_p(l)=mu_n+sqrt(var_n)*randn;
end
%x_est=mean(x_p);

mu(n)=mu_n;

end
tEnd = toc(tStart)