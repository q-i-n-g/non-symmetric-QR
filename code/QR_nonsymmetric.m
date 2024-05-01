%求出非对称实矩阵进行实schur分解Q'TQ,并求出其特征值
function [i,A,Q,L] = QR_nonsymmetric(A)
n = length(A);
L=[];%特征值数组，当前为空
if(n==1)%若矩阵为1阶，直接返回
    return
end
[A,P] = hesseborg(A);%对矩阵A进行上hessenborg化
Q = P;%得到正交矩阵
m = n;
for i = 1:10000%最大迭代次数限制
    [A,m,l] = judge(A,m);%判断收敛性，并将已经收敛的置零
    if(m<=1)%m<=1时，整个矩阵为拟上三角矩阵，完成，退出
        return
    else
        if(m-l+1>2)%若待处理的不可约的hessenborg阵维数大于2，采用双重步QR迭代
            [P,A(l:m,l:m)]=double_QR(A(l:m,l:m));
            Q(l:m,:) = P*Q(l:m,:);%更新Q
            %更新A
            if(l>1)
                A(1:l-1,l:m)=A(1:l-1,l:m)*P';
            end
            if(m+1<=n)
                A(l:m,m+1:n) = P*A(l:m,m+1:n);
            end
        else%若待处理的不可约的hessenborg阵维数等于2，特殊处理
           [A(l:m,l:m),m,P,L]=one_QR(A(l:m,l:m),m,L);
           if(P == eye(2))%若正交矩阵为单位阵，不用更新(特征值为复数的情况)
               ;
           else
               Q(l:m,:) = P*Q(l:m,:);%更新Q
               %更新A
               if(l>1)
                   A(1:l-1,l:m)=A(1:l-1,l:m)*P';
               end
               if(m+1<=n)
                   A(l:m,m+1:n) = P*A(l:m,m+1:n);
               end
           end
        end
    end
end
