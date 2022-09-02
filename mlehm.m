function [ll,lj] = mlehm(param,inp)

theta = param(1,1);
R = param(2,1);
beta = 0.9;

invest = inp.invest;
age = inp.age;
nj = inp.nj;
phat = inp.phat;
avec = [1:1:5]';

transmat0 = zeros(5,5);
for i = 1:4
    transmat0(i,i+1) = 1;
end
transmat0(5,5) = 1;

transmat1 = zeros(5,5);
for i = 1:5
    transmat1(i,1) = 1;
end

phatmat = repmat(phat,[1,5]);

v_term1 = inv(eye(5)-beta*(1-phatmat).*transmat0 - beta*phatmat.*transmat1);
v_term2 = (ones(5,1) - phat).*(theta.*avec + 0.5772.*(ones(5,1))-log(ones(5,1)-phat)) + phat.*(R + 0.5772*ones(5,1)-log(phat));

v = v_term1*v_term2;

p = exp(R*ones(5,1)+beta*transmat1*v)./(exp(theta*avec+beta*transmat0*v)+exp(R*ones(5,1)+beta*transmat1*v));


lj = zeros(nj,1);

for j = 1:nj
    aj = age(j,1);
    ij = invest(j,1);
    
    lj(j,1) = ij*(log((p(aj,1)))+(1-ij)*log((1-p(aj,1))));
end

ll = -sum(lj);
end
