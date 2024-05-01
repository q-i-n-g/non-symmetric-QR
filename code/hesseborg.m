%将方阵A正交相似变化为上hessenborg阵
function [Hes,Q] = hesseborg(A)
n = length(A);
Q = eye(n);
I = eye(n);
Hes = A;
%对第1列到第n-2列依次做househoulder变化，并做对应行变化
for i = 1:n-2
    x = Hes(i+1:n,i);
    [v,b] = house(x);
    k = length(v);
    %H=I-b*v*vT，Q=HQ，将H拆开与Q相乘，可以更快速(O(n^3)到O(n^2))
    Q(n-k+1:n,:) =Q(n-k+1:n,:)-b*v'*((v)*Q(n-k+1:n,:));%累计正交矩阵Q
    %更新hes，最终得到上hessenborg阵
    Hes(n-k+1:n,:) = Hes(n-k+1:n,:) - b*v'*((v)*Hes(n-k+1:n,:));
    Hes(:,n-k+1:n) = Hes(:,n-k+1:n) - b*(Hes(:,n-k+1:n)*v')*v;
end

