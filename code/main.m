A = randi([1,1000],50,50);
n = length(A);
B = A;
[time,A,Q,L] = QR_nonsymmetric(A);%调用获得A的实schur分解，L中包含了A的复数特征值
k = length(L);
i=1;
while(i<=n)%遍历schur分解的对角元，得到A的实数特征值
    if(i==n)
        if(A(i,i-1)==0)
            k = k+1;
            L(k) = A(i,i);
        end
        break;
    else
        if(A(i+1,i)==0)
            k = k+1;
            L(k) = A(i,i);
            i = i+1;
        else%跳过次对角元不为0的（共轭复特征值，已经加入L中）
            i = i+2;
        end
    end
end
%误差，残差计算
errorE_2norm = norm(Q*B*Q'-A)/norm(B);
errorF_2norm = norm(Q*Q'-eye(n));
errorvalue_2norm = norm(sort(L','descend')-sort(eig(B),'descend'),'inf');
disp(errorE_2norm);
disp(errorF_2norm);
disp(errorvalue_2norm);

