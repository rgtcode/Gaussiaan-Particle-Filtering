function [mu, tEnd] = SIR_UNGM(x,y)

%x = 0.1; % initial actual state
x_N = 1; % Noise covariance in the system (i.e. process noise in the state update, here, we'll use a gaussian.)
x_R = 1; % Noise covariance in the measurement (i.e. the Quail creates complex illusions in its trail!)
T = 100; % duration the chase (i.e. number of iterations).
N = 10; % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.

V = 2; %define the variance of the initial esimate

x_P = []; % define the vector of particles

% make the randomly generated particles from the initial prior gaussian distribution
tStart = tic;
for i = 1:N
    x_P(i) = x(1) + sqrt(V) * randn;
end

mu = zeros(1,N);
for t = 1:T
   % x(t+1) = 0.5*x(t) + 25*x(t)/(1 + x(t)^2) + 8*cos(1.2*(t-1)) +  sqrt(x_N)*randn;
   % y = x^2/20 + sqrt(x_R)*randn;
    %Here, we do the particle filter
    for i = 1:N
        
        x_P_update(i) = 0.5*x_P(i) + 25*x_P(i)/(1 + x_P(i)^2) + 8*cos(1.2*(t-1)) + sqrt(x_N)*randn;
        %with these new updated particle locations, update the observations
        %for each of these particles.
        y_update(i) = x_P_update(i)^2/20;
      
        P_w(i) = (1/sqrt(2*pi*x_R)) * exp(-(y(t+1) - y_update(i))^2/(2*x_R));
    end
    
    % Normalize to form a probability distribution (i.e. sum to 1).
    P_w = P_w./sum(P_w);
    
    %% Resampling: From this new distribution, now we randomly sample from it to generate our new estimate particles
  
    for i = 1 : N
        x_P(i) = x_P_update(find(rand <= cumsum(P_w),1));
    end
    
    %The final estimate is some metric of these final resampling, such as
    %the mean value or variance
    x_est = mean(x_P);

    % Save data in arrays for later plotting
    %x_out = [x_out x];
    %y_out = [y_out y];
    %x_est_out = [x_est_out x_est];
    mu(t) = x_est;
end
tEnd = toc(tStart)