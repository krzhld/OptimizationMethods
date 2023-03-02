H=[10 20 30 40 50 60 70 80 90 100];
X=ones(15,1);
for i=[10 50 100 500 1000 5000 10000 50000 100000 500000]
    D=diag([1,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ,i]);
    [Q,R]=qr(rand(15));
    A=Q*D*Q';
    b=A*X;
    writematrix(A,'C:\Users\Acer\Documents\6 семестр\Метопты\1\compare_methods\A.txt','WriteMode','append','Delimiter','tab');
    writematrix(b,'C:\Users\Acer\Documents\6 семестр\Метопты\1\compare_methods\b.txt','WriteMode','append','Delimiter','tab');
end

format long;
resultRotation=fopen('C:\Users\Acer\Documents\6 семестр\Метопты\1\compare_methods\compare methods\compare methods\rotationresult.txt','r');
resultGauss=fopen('C:\Users\Acer\Documents\6 семестр\Метопты\1\compare_methods\compare methods\compare methods\gaussresult.txt','r');
formatSpec = '%f';
allXR=fscanf(resultRotation,formatSpec);
allXG=fscanf(resultGauss,formatSpec);
fclose(resultRotation);
fclose(resultGauss);

normX=norm(X);

N=15;
colsr = floor(length(allXR)/N);
colsg = floor(length(allXG)/N);
XXXrot = reshape(allXR(1:N*colsr), [N,colsr]);
XXXgau = reshape(allXG(1:N*colsg), [N,colsg]);

H=[10 50 100 500 1000 5000 10000 50000 100000 500000];

for i=1:10
    rel_error_rot(i)=abs(norm(XXXrot(:,i)-X)/norm(X));
    rel_error_gau(i)=abs(norm(XXXgau(:,i)-X)/norm(X));
end

loglog(H,rel_error_rot);
hold on;
loglog(H,rel_error_gau);
title('Зависимость относительной погрешности от числа обусловленности');
xlabel('Число обусловленности');
ylabel('Относительная погрешность');

legend({'Rotation method','Gauss method'});