function [v0,v1] = vfi(theta1,R)

amin = 1;
amax = 5;
step = 1;
a = amin:step:amax;
a = transpose(a);
n=length(a);

beta = 0.9;

v0_0=zeros(n, 1);
v1_0=zeros(n, 1);
v0_new = zeros(n, 1);
v1_new = zeros(n, 1);
v0 = v0_0;
v1 = v1_0;
err = 100;
tol = 1e-05;
v0a1_0 = zeros(n, 1);
v1a1_0 = zeros(n, 1);
v0a1_1 = zeros(n, 1);
v1a1_1 = zeros(n, 1);

while err>tol
    v0a1_0(1:end-1) = v0(2:end);
    v0a1_0(end) = v0(end);
    v1a1_0(1:end-1) = v1(2:end);
    v1a1_0(end) = v1(end);
    
    v0a1_1 = transpose([v0(1), v0(1), v0(1), v0(1), v0(1)]);
    v1a1_1 = transpose([v1(1), v1(1), v1(1), v1(1), v1(1)]);
    
    
    v0_new = theta1*a + beta*(0.5772+log(exp(v0a1_0) + exp(v1a1_0)));
    v1_new = R + beta*(0.5772+log(exp(v0a1_1) + exp(v1a1_1)));
    err = max(max(abs(v0_new-v0)), max(abs(v1_new-v1)));
    fprintf('Error is:  %6.2e \n',err)
    v0=v0_new;
    v1=v1_new;
end