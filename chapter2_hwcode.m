%% dichotomy example root(x+lnx) epslion=5e-4 1
clc,clear
a=1/2;b=1;epslion = 5e-4;f=@(x)x+log(x);
x=a:0.01:0.75;y=f(x);plot(x,y)
title('bindiv')
k=0;errork = [];
hold on
u = fsolve(f,1);
while b-a>=2*epslion
    k=k+1;
    mid = (b+a)/2;
    errork = [errork,abs(mid-u)];
    t = f(mid);
    if t*f(a)<0
        b = mid;
    elseif t*f(b)<0
        a = mid;
    else
        break
    end
    plot(x,y)
    set(gca,'XAxisLocation','origin') 
    set(gca,'YAxisLocation','origin')
    text(mid,t,'*')
    fill(mid,t,'r')
    pause(1)
    clf
    hold on
    % text(mid,t,num2str(k),'Color','red')
end
close
disp((a+b)/2)
disp((a+b)/2-u);
plot(1:k,log10(errork))
title('误差曲线图log(err)')
%% dichotomy example root(e^x+x-3) epslion=5e-4 2
clc,clear
a=1/2;b=1;epslion = 5e-4; f2 = @(x)exp(x)+x-3;
errork=[];u2 = fsolve(f2,1);k=0;
while b-a>=2*epslion
    mid = (b+a)/2;
    k=k+1;
    errork = [errork,abs(mid-u2)];
    t = f2(mid);
    if t*f2(a)<0
        b = mid;
    elseif t*f2(b)<0
        a = mid;
    else
        break;
    end
end
disp((a+b)/2)
err_ter=(a+b)/2-u2;
disp(err_ter);
plot(1:k+1,[log10(errork),log10(err_ter)])
title('误差曲线图log(err)')
%% stationary point recur example 3
clc,clear
x0 = 1/2; f3=@(x)exp(-x);x1=f3(x0); epslion = 5e-4;
x = 0.4:0.001:0.7;
k=1;f3_t=@(x)x+log(x);
u3 = fzero(f3_t,x0);
err_k = abs(x1-u3);
lambda = @(x)1/(1+exp(-x));
while abs(x1-x0)>epslion
    k=k+1;
    x0 = x1;
    %x1 = f3(x0);
    x1 = (1-lambda(x0))*x0+lambda(x0)*f3(x0);
    err_k = [err_k,abs(x1-u3)];
    axis([0.5 0.6 -0.2 0.2])
    plot(x,f3_t(x))
    set(gca,'XAxisLocation','origin') 
    set(gca,'YAxisLocation','origin')
    text(x1,f3_t(x1),'*')
    fill(x1,f3_t(x1),'r')
    pause(1)
    clf
    hold on
end
close
figure(2)
plot(1:k,log10(err_k))
title('log误差曲线图')
disp(x1)
%% stationary point recur example 4
clc,clear
x0 = 1/2; f4=@(x)log(5-x);x1=f4(x0); epslion = 5e-4;
g = @(x)exp(x)+x-5;
x = 1:0.001:2; plot(x,f4(x)) 
hold on
plot(x,x,'r')
legend('y=ln(5-x)','y=x')
while abs(x1-x0)>epslion
    x0 = x1;
    x1 = f4(x0);
end
disp(x1)
sol = fzero(g,x0);
disp(x1-sol)
disp(x1-x0)
%% solve f(x)=sqrt(x),x=115 f'(x)=1/2sqrt(x) example 5
% 他是谁的根
x0 = 11; f5 = @(x) x^2-115;x1 = x0-f5(x0)/(2*x0); epslion = 5e-5;
k=1;err_k = abs(x1-sqrt(115));
while abs(x1-x0)>epslion
    k=k+1;  
    x0 = x1;
    x1 = x0-f5(x0)/(2*x0);
    err_k = [err_k, abs(x1-sqrt(115))];
end
disp(x1)
disp(x1-sqrt(115))
plot(1:k,log10(err_k))

title('log误差曲线图')
disp(x1)
%% iter -1/6<c<0, 6
clc,clear
c = -2/(3*7^(2/3))+0.0001;%c=-1/6
x0 = 2.1; f6=@(x)x+c*(x.^3-7);x1=f6(x0); epslion = 5e-4;
% x = 1.8:0.001:2; plot(x,f6(x)) 
% hold on
% plot(x,x,'r')
% legend('y=x-x^3-7/6','y=x')
k=1;
while abs(x1-x0)>epslion && k<300
    x0 = x1;
    x1 = f6(x0);
    k=k+1;
end
disp(x1)
disp(['与真解的误差是',num2str(x1-7^(1/3))])
%% suitable recur newton method 7
clc,clear
% x0 = 0.5; f7=@(x)1/(x+1); epslion = 5e-4;
% u7 = fsolve(@(x)x^2+x-1,1/2);x1 = f7(x0);
% err_k=abs(x1-u7);k=1;
% while abs(x1-x0)>epslion
%     x0 = x1;
%     x1 =f7(x0);
%     k=k+1;
%     err_k = [err_k, abs(x1-u7)];
% end

% 
x0 = 0.001; f7=@(x)x^2+x-1; x1 = x0-f7(x0)/(2*x0+1); epslion = 5e-4;
u7 = fsolve(f7,1/2);err_k=abs(x1-u7);k=1;
while abs(x1-x0)>epslion
    x0 = x1;
    x1 = x0-f7(x0)/(2*x0+1);
    k=k+1;
    err_k = [err_k, abs(x1-u7)];
end
plot(1:k,log10(err_k))
title('log误差曲线图')
disp(x1)
disp(x1-u7)

%% newton method example 8
clc,clear
x0 = 0.409; f8=@(x)4*x.^3-8*x.^2+5*x-1; x1 = x0-f8(x0)/(12*x0^2-16*x0+5); epslion = 5e-8;
x=-4:0.1:5;k=1;
u8 = fsolve(f8,x0);err_k=abs(x1-u8);
plot(x,f8(x))
set(gca,'XAxisLocation','origin') 
set(gca,'YAxisLocation','origin')
while abs(x1-x0)>epslion
    k=k+1;
    x0 = x1;
    x1 = x0-f8(x0)/(12*x0^2-16*x0+5);
    err_k = [err_k, abs(x1-u8)];
end
disp(x1)
plot(1:k,log10(err_k))
title('log误差曲线图')
disp(x1-x0)
k
%% newton secant method 9
clc,clear
x0 = 1;x1=3/2; f8=@(x)x^3+x-3; x2 = x1-f8(x1)/(f8(x1)-f8(x0))*(x1-x0); epslion = 5e-4;
u8 = fsolve(f8,x0);k=1;
while abs(x1-x2)>epslion
    x0 = x1;
    x1 = x2;
    x2 = x1-f8(x1)/(f8(x1)-f8(x0))*(x1-x0);
    k=k+1;
end
disp(x2)
disp(x2-u8)