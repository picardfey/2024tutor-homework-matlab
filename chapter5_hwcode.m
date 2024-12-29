%% 2


%% 5 Gauss integer
clc,clear
f = @(x)1/(1+x^4);
syms x
%two points
% 换元 
ltx = @(t)1.5*t-0.5;
ans1 = (f(ltx(-1/sqrt(3)))+f(ltx(1/sqrt(3))))*1.5;
ans2 = (5*f(ltx(-sqrt(3/5)))+8*f(ltx(0))+5*f(ltx(sqrt(3/5))))/9*1.5;
symf = symfun(f(x), x);
u = double(int(symf, -1, 2));
%% 6
clc,clear
time = 5*1:10;
velo = [20 16 12 8 11 14 17 21 16 11];
dis = 5*(sum(velo)-(velo(1)+velo(end))/2);
%向前差商/向后差商
a1 = (velo(2:end)-velo(1:end-1))/5;
