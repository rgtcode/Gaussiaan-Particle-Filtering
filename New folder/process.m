function out=process(x,i)
alpha=0.5;
beta=25;gamma=8;

%x_N=1;  %Noise variance at measurement update

out=alpha*x+beta*(x/(1+x^2))+gamma*cos(1.2*(i-1));
end