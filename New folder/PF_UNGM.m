%% Initialization
x_initial=0.1; %initial state
x_N=1;  %Noise variance at measurement update
x_R=1;  %Noise variance at time update
N=100;  %No. of states
M=100;   %No. of particles
%x_p_update=zeros(1,M);
%y_update=zeros(1,M);
alpha=0.5;
beta=25;gamma=8;
v=1;%variance around the initial state

x = zeros(1,N+1);
x(1) = x_initial;
y = zeros(1,N+1);
y(1) = x(1)^2/20+sqrt(x_R)*randn;

for i=2:N+1 
    x(i)=alpha*x(i-1) + beta*(x(i-1)/(1+x(i-1)^2)) + gamma*cos(1.2*(i-2)) + sqrt(x_N)*randn;
    y(i)=x(i)^2/20+sqrt(x_R)*randn;
end

%% Calling the Filters
[mu_gpf, time_gpf] = GPF_UNGM(x,y);
[mu_sir, time_sir] = SIR_UNGM(x,y);
[mu_ekf, time_ekf] = EKF_UNGM(x,y);

MSE_gpf = sum((mu_gpf- x(2:end)).^2,1);
MSE_sir = sum((mu_sir- x(2:end)).^2,1);
MSE_ekf = sum((mu_ekf- x(2:end)).^2,1);

%% MSE plot
figure(1)
plot(1:N, MSE_gpf,'-o',1:N,MSE_sir,'-+',1:N,MSE_ekf,'-d');
legend('GPF','SIR','EKF')
xlabel('realizations')
ylabel('MSE')

%% Computation Time

comp_time= [time_sir, time_gpf, time_ekf];
names = ['SIR', 'GPF', 'EKF'];
figure(2)
h = bar(comp_time);
l = cell(1,3);
l{1}='SIR'; l{2}='GPF'; l{3}='ekf';  
set(gca,'xticklabel', l) 
ylabel('Computation time')


