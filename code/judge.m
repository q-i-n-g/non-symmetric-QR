%收敛性判定函数
%从i行开始向上走，获得最大非负整数m，m行以下的矩阵为拟上三角阵
%获得最小非负整数l，使得l行到m行，l列到m列的子矩阵为不可约的hessenborg阵
function [A,m,l] = judge(A,i)
m = i;
l = i;
while(l>1)
    if ( abs(A(l,l-1)) <= eps*(abs(A(l,l))+abs(A(l-1,l-1))) )
        A(l,l-1) = 0;%将次对角线hi-1,i上绝对值<=u*(|hi,i|+|hi-1,i-1|)的元素赋为零
        if(m == l) %m=l，说明此时还未找到最大拟上三角阵，继续向上找           
            m = m-1;
            l = l-1;            
        else%m!=l，但当前次对角元为0，说明此处hessenborg子矩阵开始不可约，停止
            break;
        end
    else%次对角线元素不满足收敛条件时，说明找到最大拟上三角阵，m不再变化，l继续向上
            l = l-1;
    end
end
