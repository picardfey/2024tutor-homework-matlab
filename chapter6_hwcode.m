%% 1
% 欧拉法
clc,clear
y0 = 1;
N = 10; h =2/N;f = @(x, y)x^2-y;
y_1 = zeros(1, N+1); y_1(1)=y0;
y_modify = y_1;
x = (0:N)*h;  %避免浮点数误差的累积 0:h:2
for i=1:N
    fi = f(x(i), y_1(i));
    fi_modify =  f(x(i), y_modify(i));  %分开写
    y_pred = y_modify(i)+h*fi_modify;
    y_1(i+1) = y_1(i)+h*fi;
    y_modify(i+1) = y_modify(i)+h/2*(fi_modify+f(x(i+1),y_pred));
end
[xs,ys] = ode45(f, [0 2], 1);
plot(x, y_1,'b--', x,y_modify,'r-', xs, ys)
legend('欧拉法','改进欧拉法','真解')
% 分割数据
half_idx = ceil(N/ 2+1); 

% 打印结果表格
fprintf('x_i(左) | 欧拉法 (左)   |改进欧拉法 (左) ||  x_i (右)| 欧拉法 (右)    | 改进欧拉法 (右)\n');
fprintf('----------------------------------------------------------------------------------------\n');
for i = 1:half_idx
    if i + half_idx <= N+1
        fprintf('%6.4f  |  %10.6f  |  %10.6f  ||  %6.4f  |  %10.6f  |  %10.6f\n', ...
            x(i), y_1(i), y_modify(i), ...
            x(i + half_idx), y_1(i + half_idx), y_modify(i + half_idx));
    else
        fprintf('%6.4f  |  %10.6f  |  %10.6f\n', ...
            x(i), y_1(i), y_modify(i));
    end
end
%% 2
clc,clear
N = 10;h = (2-1)/N;
x = 1+(0:N)*h;
y_1 = zeros(1, N+1);
y_1(1) = 0;y_trap = y_1;
yx = @(x)x.^2.*(exp(x)-exp(1));
f = @(x, y)2*y/x+x^2*exp(x);
for i=1:N
    fi = f(x(i), y_1(i));
    y_1(i+1) = y_1(i)+h*fi;
    % 隐式方法，解非线性方程 第二章
    func = @(u)u-y_trap(i)-h/2*(f(x(i), y_trap(i))+f(x(i+1), u));
    y_trap(i+1) = fzero(func, y_1(i+1));
end
y_true = yx(x);
%y_true = ode45(f, [1 2], 0);
err_1 = y_true - y_1; % 欧拉法误差
err_trap = y_true - y_trap; % 梯形法误差
%打印表格
fprintf('  x_i     | 欧拉法    | 梯形法    | 精确解    | 欧拉法误差  | 梯形法误差\n');
fprintf('--------------------------------------------------------------------------\n');
for i = 1:N+1
    fprintf('%6.4f  | %10.6f | %10.6f | %10.6f | %10.6f | %10.6f\n', ...
            x(i), y_1(i), y_trap(i), y_true(i), err_1(i), err_trap(i));
end

% 画图
figure(1)
plot(x, y_1, 'b--', x, y_trap, 'r-', x, y_true, 'k-');
legend('欧拉法', '改进欧拉法', '解析解');
title('欧拉法 vs 改进欧拉法 vs 解析解');
xlabel('x'); ylabel('y');

figure(2)
plot(x, err_1, 'b--', x, err_trap, 'r-');
legend('欧拉法误差', '改进欧拉法误差');
title('误差对比');
xlabel('x'); ylabel('误差');

%%  3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear
for j=1:2
    N = 15*2^(j-1); h = 3/N;
    x = 0+(0:N)*h;y = zeros(1, N+1);
    y(1) = 0.5;
    f = @(x,y)-9*y;
    for i=1:N
        fi = f(x(i), y(i));
        y_pred = y(i)+h* fi;
        y(i+1) = y(i)+h/2*(f(x(i+1),y_pred)+fi);
    end
    hold on
    plot(x, y, '-r')
    y_true = 0.5*exp(-6*x);
    err = norm(y_true(1:2^j:end)-y(1:2^j:end), 2);
    disp(['步长为0.', num2str(3-j), '时，整体误差为',num2str(err)])
end
%%   4
clc, clear
f = @(x, y)-x*y+2*x^2+2*x;
N = 10; h =1/N;
x = 0+(0:N)*h;  %避免浮点数误差的累积
y = zeros(1, N+1); y(1)=1;
y_gill = y;
for i=1:N
    k1 = f(x(i), y(i));k1_gill = f(x(i), y_gill(i));
    k2 = f(x(i)+h/2, y(i)+h*k1/2);k2_gill = f(x(i)+h/2, y_gill(i)+h*k1_gill/2);
    % k1, k2取法相同
    k3 = f(x(i)+h/2, y(i)+h*k2/2);k3_gill = f(x(i)+h/2, y_gill(i)+(sqrt(2)-1)/2*h*k1_gill+(2-sqrt(2))/2*h*k2_gill);
    % gill公式较于RK4斜率的取法不同
  
    k4 = f(x(i+1), y(i)+h*k3);k4_gill = f(x(i+1), y_gill(i)-sqrt(2)/2*h*k2_gill+(2+sqrt(2))/2*h*k3_gill);
    y(i+1) = y(i)+h*(k1+2*k2+2*k3+k4)/6;
    y_gill(i+1) = y_gill(i)+h*(k1_gill+(2-sqrt(2))*k2_gill+(2+sqrt(2))*k3_gill+k4_gill)/6;
end
plot(x, y, '-r',x, y_gill, '*b')
legend('RK4','Gill')
u_true = ode45(f, [0 1], 1).y;
%% 5
clc,clear
N = 20;h = 2/N;
x = 0+h*(0:N);y = zeros(1, N+1);
y(1)=1;
f = @(x, y)-y+x;
for i=1:N
    fi = f(x(i), y(i));
    y(i+1) = y(i)+h*fi;
end
plot(x, y)
legend('数值解曲线')
x_vertify = [0.5, 1.0, 1.5];
y_true = @(x)x-1+2*exp(-x);
index = [1+0.5/h, 1+1.0/h, 1+1.5/h];
% 局部截断误差
taun = abs(y_true(x_vertify)-y_true(x_vertify-h)-h*(f(x(index-1), y_true(x_vertify-h))));
disp(['局部截断误差分别是', num2str(taun)])
% 整体截断误差
e_n = abs(y_true(x_vertify)-y(index));
disp(['整体截断误差分别是',num2str(e_n)])
%% 6
N = 20; h=1/N;
x = 0+h*(0:N);y = zeros(1, N+1);
y(1) = 1;
f = @(x, y)-y;
for i=1:N
    fi = f(x(i), y(i));
    y_pred = y(i)+h*fi;
    y(i+1) = y(i)+h/2*(fi+f(x(i+1),y_pred));
end
plot(x, y)
% 取蓟县,exp(-nh)


