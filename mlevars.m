function vars_matrix = mlevars(param0,stpsize,inp)

invest = inp.invest;
age = inp.age;
nj = inp.nj;

pi = 3.141592;

llobji = zeros(nj,1);
llobji0 = zeros(nj,1);
deltalni = zeros(nj,1);
deltalni2 = zeros(nj,2);

for coeff = 1:2
    theta10 = param0(1,1);
    R0 = param0(2,1);
    param_org = [theta10;R0];
    
    param = param0;
    param(coeff,1) = param0(coeff,1)*(1+stpsize); % Change param values
    
    theta1 = param(1,1);
    R = param(2,1);
    param_new = [theta1;R];
    
    
    [~,llobji(:,1)] = mle(param_new,inp);
    
    [~,llobji0(:,1)] = mle(param_org,inp);
    
    
    deltalni(:,1) = (llobji(:,1) - llobji0(:,1))./(stpsize*param0(coeff,1));
    deltalni2(:,coeff) = deltalni;
end

stack_indiv_var = zeros(2,2,nj);

for i = 1:nj
    stack_indiv_var(:,:,i) = deltalni2(i,:)'*deltalni2(i,:);
end

sum_indiv_var = sum(stack_indiv_var,3);

vars_matrix = inv(sum_indiv_var);
