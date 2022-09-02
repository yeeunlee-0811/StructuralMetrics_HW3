clear

%% Q3: valFun iter.

theta = -1;
R = -3;
beta = 0.9;

[v0,v1] = vfi(theta,R);

%% Q5: 
load data.asc

%%
age = data(:,1);
invest = data(:,2);
nj = size(age,1);

inp.invest = invest;
inp.age = age;
inp.nj = nj;

%%
llobj = @(param)mle(param,inp);

param0 = [-1;-3];

options = optimset('Display','iter','PlotFcns',@optimplotfval);

[mleparam,logval] = fminsearch(llobj,param0,options);

%% 
stp = 1e-5;
mlevarmat = mlevars(mleparam,stp,inp);

mlese = sqrt(diag(mlevarmat));

%% 
phat = zeros(5,1);

for i = 1:5
    iage = age ==i;
    ageinv = invest(iage);
    iiage = ageinv ==1;
    phat(i,1) = sum(iiage)./sum(iage);
end

inp.phat = phat;

%%
hmllobj = @(param)mlehm(param,inp);

param0 = [-1;-3];

options = optimset('Display','iter','PlotFcns',@optimplotfval);

[mleparam_hm,logvalhm] = fminsearch(hmllobj,param0,options);

