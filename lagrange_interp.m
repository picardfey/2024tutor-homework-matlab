function L = lagrange_interp(x, y)
    % x: 插值节点的横坐标数组
    % fx: 插值节点的函数值数组
    % x_eval: 要求的插值点
    % L: x_eval 处的插值结果

   n = length(x);
   syms X
   L = 0;
   for i=1:n
       Li = 1;
       for j=1:n
            if j~=i
                Li = Li*(X-x(j))/(x(i)-x(j));
            end
       end
       L = L+y(i)*Li;
   end