%% horner algorithm example 5
% 点击运行节，只运行当前节块
f = 100*ones(1,4);
x = [0.2 0.4 0.7 0.9]; %[0.2 0.4,0.7,0.9] 
for i = 100:-1:1
    f = f.*x+i-1;
end
disp(f+1)
disp(1+x./(1-x).^2)
%% function approximate sheme who best 6
x = 1.400:0.0001:1.4142;
y = (x-1).^6;
y1 = abs(1./(x+1).^6-y);
y2 =abs(1./(3+2*x).^3-y);
y3 = abs((3-2*x).^3-y);
% y4 = 99-70*x;
y5 = abs(1./(99+70*x)-y);
plot(x,y1,'-b',x,y2,'*r',x,y3,'dg',x,y5,'.k')
legend('y1','y2','y3','y5')
%% your computer disposition huge+small num 7
% S=0;n=1;
% while S+1/n>S
%     S = S+1/n;
%     n=n+1;
% end
% k=0;sub=0;
% for i=1:4
%     k=k+1/(n-i);
%     sub = [sub,k];
% end
% disp(S-sub)
% 太折磨人了，其实完全可以这么做计算数量级多大会被吃
% solve k:2^k<=ln(2^(53-k))+0.5771<=2^(k+1)
% hence n=2^(53-k)=2^48
% 16个数量级,你的电脑要跑到10^14去
%% stable recur algorithm int(x^n*e^x,0,1) 8
n = 0:6;
f = @(x) x.^n.*exp(x);
z = integral(f,0,1,'ArrayValued',true);
disp(z)

Intk = (exp(1)+1)/200;
res = zeros(1,7);
for n=1:99
    Intk = (exp(1)-Intk)/(100-n);
    if n>92
        res(n-92)=Intk;
    end
end
disp(flip(res))
%% forward euler solve 1st constant ode 9
e0 = 0.01;h=0.1:0.1:0.3;
for i=1:100
    e0 = e0-6*h.*e0;
end
disp(e0)

n = 30;hs=0.5;
x = 0:hs:n*hs;
y=zeros(1,n+1);y(1)=1;
for i=1:n
    y(i+1) = y(i)-6*hs*y(i);
end
plot(x,exp(-6*x),'-r',x,y,'*b')
legend('truth','proximate')
