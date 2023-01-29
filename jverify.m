function flag = jverify(u)

global N
global t
r = 1e-3*rand(N,1);

[f0,j0] = calFJ(u);
figure(1)
[f1,~] = calFJ(u + r);
feval0 = f1 - f0 - j0*r;
for iter = 2:10
    ur = u + iter*r;
    [f_i,~] = calFJ(ur);
    feval_i = f_i - f0 - j0*(iter*r);
    feval_i = sqrt(feval_i./feval0);
    plot(t,feval_i);
    hold on
end

flag = 1;
end