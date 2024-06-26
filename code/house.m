%求将向量x镜像变换到e1方向([1 0 0 ...]T)的正交矩阵H = I - b*v*vT
function [v,b] = house(x)
n = length(x);
y = norm(x,"inf");
x = x/y;%避免计算分量平方时上溢
z = x(2:n)'*x(2:n);
v(2:n)=x(2:n);
if (z==0)%已经是在e1方向/发生下溢
    b = 0;
else
    a = (x(1)^2+z)^0.5;%a为x的二范数
    if(x(1)<=0)%x(1)<0时，取v=x1-||x||2是避免了减法相消的好的选取
        v(1)=x(1)-a;
    else%x(1)>0时，取v=x1-||x||2存在减法相消，为了避免减法相消，拆开计算v
        v(1) = -z/(x(1)+a);
    end
    b = (2*v(1)^2)/(z+v(1)^2);%b=2/(vT*v)，将x拉回原本长度
    v = v/v(1);
end

