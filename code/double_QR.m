%双重步位移的QR迭代
function [Q,A]=double_QR(A)
n = length(A);
m = n-1;
%两次位移后M=H^2-sH+tI,s=hmm+hnn(共轭复根之和),t=det(H(n-1:n,n-1:n)(共轭复根之积))
s = A(m,m)+A(n,n);
t = A(m,m)*A(n,n)-A(m,n)*A(n,m);
x = A(1,1)*A(1,1)+A(1,2)*A(2,1)-s*A(1,1)+t;
y = A(2,1)*(A(1,1)+A(2,2)-s);
z = A(2,1)*A(3,2);
Q = eye(n);
I = eye(3);
for k= 0:n-3%bulge—chasing过程
    [v,b] = house([x,y,z]');
    H = I - b*(v'*(v));%Q由q1唯一确定，因此求出q1
    Q(k+1:k+3,:) = H*Q(k+1:k+3,:);
    %bluge—chasing过程，将矩阵重新化为上hessenborg阵
    q = max([1,k]);
    A(k+1:k+3,q:n) = H*A(k+1:k+3,q:n);
    r = min([k+4,n]);
    A(1:r,k+1:k+3) = A(1:r,k+1:k+3)*H';
    x = A(k+2,k+1);
    y = A(k+3,k+1);
    if(k<n-3)
        z = A(k+4,k+1);
    end
end
%bluge—chasing的最后一步
[v,b] = house([x,y]');
H = eye(2) - b*(v'*(v));
Q(n-1:n,:) = H*Q(n-1:n,:);
A(n-1:n,n-2:n) = H*A(n-1:n,n-2:n);
A(1:n,n-1:n) = A(1:n,n-1:n)*H;
