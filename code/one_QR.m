%处理双重步迭代剩下的两阶矩阵
function [A,m,Q,L] = one_QR(A,m,L)
a = 1;
b = -(A(1,1)+A(2,2));
c = A(1,1)*A(2,2)-A(1,2)*A(2,1);
lamda1 = (-b+(b^2-4*a*c)^0.5)/(2*a);%计算这个两阶矩阵的特征值
lamda2 = (-b-(b^2-4*a*c)^0.5)/(2*a);
k = length(L);
Q = eye(2);
if(imag(lamda1)==0)%如果特征值为实数，则进行带位移的QR
    lamda = A(2,2);
    A = A - lamda*eye(2);
    x = A(:,1);
    [v,b] = house(x);
    H = eye(2) - b*v'*v;
    A = H*A*H'+lamda*eye(2);
    Q = H;
else%如果特征值为复数，则不对这个二阶矩阵进行处理
    L(k+1) = lamda1;%将复数特征值放入特征值数组中（方便之后求特征值）
    L(k+2) = lamda2;
    m = m-2;
end
