function PF_BOT

x_initial =[-0.05 0.001 0.7 -0.055];                                               %initial state
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

x = zeros(N,4);
y = zeros(1,25);
x(1,:) = x_initial;
y(1) = normrnd(atan(x(1,3)/x(1,1)),x_R);
for i = 2:25
    x(i,:)= mvnrnd( phi * x(i-1,:)', gamma * x_N *(gamma') );                          %Process equation
    y(i)= normrnd( atan(x(i,3)/x(i,1)), x_R);                                     %Observation equation
end

for runs= 1:1
%% GPF
mu_gpf = GPF_BOT(x,y);
%% SIR
mu_sir = SIR_BOT(x,y);
%% EKF
mu_ekf = EKF_BOT(x,y);

%% Plotting
figure(runs)
scatter(x(:,1),x(:,3), 'x','LineWidth',1.5)
hold on
scatter(mu_gpf(:,1),mu_gpf(:,3),'o','LineWidth',1.5)
hold on
scatter(mu_sir(:,1),mu_sir(:,3),'+','LineWidth',1.5)
hold on
scatter(mu_ekf(:,1),mu_ekf(:,3),'d','LineWidth',1.5)
legend('True state','GPF','SIR','EKF')
end

%% MSE plot
%figure(2)
MSE_gpf = ((mu_gpf-x(2:end,:)).^2);