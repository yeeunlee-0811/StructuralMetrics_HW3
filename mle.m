function [ll,lj] = mle(param,inp)

theta = param(1,1);
R = param(2,1);

[v0,v1] = vfi(theta,R);

rprob = exp(v1)./(exp(v0)+exp(v1));
nrprob = ones(5,1) - rprob;

invest = inp.invest;
age = inp.age;
nj = inp.nj;

lj = zeros(nj,1);

for j = 1:nj
    aj = age(j,1);
    ij = invest(j,1);
    
    lj(j,1) = ij*(rprob(aj,1))+(1-ij)*(nrprob(aj,1));
end

ll = -sum(log(lj));
end


    

