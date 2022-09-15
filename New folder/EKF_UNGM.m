function [mu, tEnd] = EKF_UNGM(x,y)

N=100;%No. of states
%x=0.1;%intial condition
p=1;%initial variance

x_R=1;
x_N=1;

x_initial=x(1);
x_est_state=[x_initial];
x_tr=x;

tStart = tic;
%Extended Kalman Filter
for i=1:N
    jacob_process=0.5+25*((1-x(i)^2)/(1+x(i)^2)^2);
    jacob_update=x(i)/10;
    
    %
    p=jacob_process*p*jacob_process'+x_N;
    S=jacob_update*p*jacob_update'+x_R;
    
    %Calculating kalman gain
    K=p*jacob_update'./S;
    
    %State prediction
    x_f=process(x(i),i);
    %Measurement prediction
    y_f=update(x_f);
    %Uodating/Correcting the predicted state using the kalman gain
    x_u=x_f+K*(y(i)-y_f);
    
    %Updating/Correcting measurement update
    y_u=update(x_u);
    p_u=(1-K*jacob_update)*jacob_process;
    p=p_u;
    %x=x_u;
    %x_est_state=[x_est_state x];
    mu(i) = x_u;
end
tEnd = toc(tStart)