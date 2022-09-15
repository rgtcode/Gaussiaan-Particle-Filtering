function mu = SIR_BOT(x,y)

                                             
x_N= (10^-6).*eye(2);                                                      %Noise variance at measurement update 
x_R= 0.005^2;                                                              %Noise variance at time update
N=24;                                                                      %No. of states
M=10000;                                                                   %No. of particles


phi =[1 1 0 0;0 1 0 0;0 0 1 1;0 0 0 1];
gamma=[0.5 0;1 0;0 0.5;0 1];

v=[0.1 0 0 0;
   0 0.005 0 0;
   0 0 0.1 0;                                                              %(3,4) element is given as 1 in the paper.
   0 0 0 0.01];                                                            %variance around the initial state

x_p = zeros(M,4);
num = 1000;
%% Samples from the prior distribution
for k = 1:num
x_pd = mvnrnd(x(1,:), v, M);  % M samples from the prior distribution
x_p = x_p + x_pd;
end
x_p = 1/num .* x_p;

%% Initialization
%y = normrnd(atan(x(3)/x(1)),x_R); %y_out

%x_out=[x];
%x_est=[x]; 
%x_est_out=[x_est];
no_of_nan = 0;
%figure(2)
p_w = zeros(1,M);
s = zeros(M,1);
%% From state 0
for l= 1:M
 s(l) =  atan(x_p(l,3)/x_p(l,1)); 
 p_w(l)= normpdf(y(1), s(l), x_R);
end
p_w = p_w ./ sum(p_w);
mu0 = p_w * x_p;
var0 = zeros(4,4);

for k=1:M
      var0= var0 + p_w(k).*((x_p(k,:)-mu0)'*(x_p(k,:)-mu0));  % variance of point estimate
end
x_p = mvnrnd(mu0,var0,M);

%% Particle Filtering
x_p1 = zeros(M,4);
%p_w  = 1/M.* ones(1,M);  %putting it at ones was giving some output.
mu = zeros(N,4);
stor_p1 = zeros(N,4);
y_p1 = zeros(M,1);
for i= 1:N
    %x= mvnrnd( phi * x', gamma * x_N *(gamma') );                          %Process equation
    %y= normrnd( atan(x(3)/x(1)), x_R);                                     %Observation equation
    prev_p_w = p_w;
    for j=1:M
        x_p1(j,:) = mvnrnd( phi * x_p(j,:)', gamma * x_N *(gamma') );
        y_p1(j) = atan(x_p1(j,3)/x_p1(j,1));
        p_w(j) = normpdf(y(i+1),y_p1(j),x_R);
    end
    
    p_w = p_w ./ sum(p_w);
    
    if isnan(p_w)
        p_w(:) = prev_p_w;
        no_of_nan = no_of_nan + 1;
    end 
%% resampling
    for l= 1:4
        x_pr = x_p1(:,l);
     for k = 1 : M
        x_p(k,l) = x_pr(find(rand <= cumsum(p_w),1));
     end
    end
    
     mu(i,:) = mean(x_p);
     
%% storing    
%x_est = mu;
%x_out=[x_out;x];
%y =[y ;y];
%x_est_out=[x_est_out; x_est];

end