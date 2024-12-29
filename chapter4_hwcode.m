x = [8, 27, 64];       % 插值节点的 x 坐标
fx = [2, 3, 4];      % 插值节点的函数值
x_eval = 40;        % 要求的插值点

L = lagrange_interp(x, fx);
% 显示结果
y_eval = double(subs(L, x_eval));
disp(['插值结果在 x = ', num2str(x_eval), ' 处的值为: ', num2str(y_eval)]);

%%
xi = [1.5 1.6 1.7];
y = [0.99749 0.99957 0.99166];
[poly,interdiff] = newton_interp(xi, y);
disp(['插值多' ...
    ['' ...
    '项式是'],char(poly)])
p = [1.55 1.65];
disp(double(subs(poly, p)))
x = 1.5:0.01:1.7;f = double(subs(poly, x));
plot(x, sin(x),x, f,'r*')
legend('sin(x)', '牛顿插值二次多项式')
%---------------------------------------------------------------------------
wucha = sin(p)-double(subs(poly, p));
disp(['误差是', num2str(wucha)])
%% 手写样条插值
x = 0:10;
y = [2.51 3.30 4.04 4.70 5.22 5.54 5.78 5.40 5.57 5.70 5.80];
A = 2*eye(11)+diag([1 1/2*ones(1,9)]',1)+diag([1/2*ones(1,9), 1]',-1);
der0 = 0.8;der10=0.2;
d = zeros(1, 9);
for i=1:9
    d(i) = 3*(y(i+2)-2*y(i+1)+y(i));
end
d = [6*(y(2)-y(1)-der0) d 6*(der10-y(11)+y(10))]';  %构建右边向量
M = A\d;  % 求出系数M
plot(x, y,'*');
spline_fun = struct();
for j=1:10   %构造线性方程组
    xlin = x(j):0.01:x(j+1);
    f = @(u)(x(j+1)-u).^3/6*M(j)+(u-x(j)).^3/6*M(j+1)+(y(j)-M(j)/6)*(x(j+1)-u)+(y(j+1)-M(j+1)/6)*(u-x(j));
    spline_fun(j).fun = f;
    hold on
    h = plot(xlin, f(xlin),'-r');
    
end
legend('插值节点', '三次样条函数')
%%
clc,clear
A1 = [1 1
    1 -2
    3 2];
b1 =[1 2 3]';
A2 = [1 -1 1
     -1 -2 2
      3 -1 -1
      1 1 3];
b2 = [2 0 1 4]';
A3 = [1 -1 1
     -1 -2 2
      3 -1 -1
      1 1 3
      2 2 2];
b3 = [2 0 1 4 3]';
x1 = A1'*A1\(A1'*b1);
%x1 = lsqminnorm(A1, b1); 
x2 = A2'*A2\(A2'*b2);
%x2 = lsqminnorm(A2, b2); 
x3 = A3'*A3\(A3'*b3);
%x3 = lsqminnorm(A3, b3); 
%%
r = [2.70 2.00 1.61 1.20 1.02]';
phi = [48 67 83 108 126]'/180*pi;
% 对式子变形
A = [ones(5, 1) r.*cos(phi)];
coeff1 = A\r;
p1 = coeff1(1);e1=coeff1(2);
err1 = norm(r-p1./(1-e1*cos(phi)),2)^2;
disp(['平方误差是',num2str(err1)])
% 或者

y = 1./r; phi = cos(phi);
x = [ones(5,1) phi];
% 最小二乘回归
coeff2 = x\y;
p2 = 1/coeff2(1);e2 = -p2*coeff2(2);


plot(phi,r, '*')
hold on
cosphi = min(phi):0.001:max(phi);
plot(cosphi, p2./(1-e2*cosphi), '-r')
% 计算平方误差
err2 = norm(r-p2./(1-e2*phi),2)^2;
disp(['平方误差是',num2str(err2)])
%%
A = [1 1 0
    0 1 1
    1 2 1
    0 0 1];
b = [2 2 3 4]';
x_hat = A\b;
norm(A*x_hat-b, 2)
%%
y = log([6 2 1 1]');
x = log([1 2 3 4]');
A = [ones(4, 1), x];
coeff = A\y;
a = exp(coeff(1))