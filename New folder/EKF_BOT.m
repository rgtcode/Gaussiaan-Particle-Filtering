function mu = EKF_BOT(x,y)
%x=[-0.05 0.001 0.7 -0.055]; %initial state
x_N= (10^-6).*eye(2);                                                      %Noise variance at measurement update 
x_R= 0.005^2;                                                              %Noise variance at time update
N=24;                                                                      %No. of states
M=1000;                                                                    %No. of particles

x_est_state=[x(1,:)];
phi =[1 1 0 0;0 1 0 0;0 0 1 1;0 0 0 1];
gamma=[0.5 0;1 0;0 0.5;0 1];
Q=10^-6*[0.5 1 0 0 ;0 0 0.5 1;0.25 0.5 0 0;0 0 0.25 0.5];

v=[0.1 0 0 0;
   0 0.005 0 0;
   0 0 0.1 0;                                                              %(3,4) element is given as 1 in the paper.
   0 0 0 0.01];                                                            %variance around the initial state

x_p_re=zeros(M, 4);
x_noisy=x(1,:);
x_tr_state=[x(1,:)];
%Actual model simultion
%{
for k=1:N
   
    x_noisy= mvnrnd( phi * x_noisy', gamma * x_N *(gamma') ); 
    y_noisy= normrnd(atan(x(3)/x(1)), x_R);
    y_tr_state(k)=y_noisy;
    x_tr_state=[x_tr_state;x_noisy];
end
%}
mu = zeros(N,4);
%Extended Kalman Filter
for i=1:N
    jacob_proces =[1 1 0 0;0 1 0 0;0 0 1 1;0 0 0 1];
    R=(x(i,1)^2+x(i,3)^2);
    jacob_update=[-x(i,3)/R 0 x(i,1)/R 0];    

     v=jacob_proces*v*jacob_proces'+Q;
     S=jacob_update*v*jacob_update'+x_R;
 
%Calculating kalman gain
     K=v*jacob_update'/S;
    
     %State prediction
      x_f =phi*x(i,:)';
      x_f=x_f';
%Measurement prediction
      y_f=normrnd(atan(x_f(3)/x_f(1)), x_R);
%Uodating/Correcting the predicted state using the kalman gain
      x_u=x_f+K'.*(y(i)-y_f);
     
%Updating/Correcting measurement update
      y_u=atan(x_u(3)/x_u(1));
      v_u=(eye(4)-K*jacob_update)*v*(eye(4)-K*jacob_update)'*(K*S*K');
      v=v_u;
      %x=x_u;
      %x_est_state=[x_est_state;x];
      mu(i,:) = x_u;
end
