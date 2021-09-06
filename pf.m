%粒子滤波

%预测方程 x(i)=sin(x(i-1))+5*x(i-1)/(x(i-1)^2+1)+Q
%观测方程 y(i)=x(i)^2+R

t=0.01:0.01:1;
x=zeros(1,100);
y=zeros(1,100);

x(1)=0.1;
y(1)=0.01^2;

for i=2:100
    x(i)=sin(x(i-1))+5*x(i-1)/(x(i-1)^2+1);
    y(i)=x(i)^2+normrnd(0,1);
end
%plot(t,x,t,y,'LineWidth',2)

n=100;
xold=zeros(1,n);
xnew=zeros(1,n);
xplus=zeros(1,n);%xplus用于存放滤波值，就是每一次后验概率的期望
w=zeros(1,n);
%设置x0(i),可以直接在正态分布采样，如果对初值有信心，也可以让所有粒子都相同
for i=1:n
    xold(i)=0.1;
    w(i)=1/n;
end

for i=2:100
    
    %预测步
    for j = 1:n
        xold(j)=sin(xold(j))+5*xold(j)/(xold(j)^2+1)+normrnd(0,0.1) %Q;
    end
    
    %更新步
    for j = 1:n
        %w(j)=w(j)*fR(y-观测)
        %正态分布
        %fR=(2*pi*R)^(-0.5)*exp(-((y(i)-xold(j)^2)^2/(2*R)))
        %w(j)=(2*pi*R)^(-0.5)*exp(-((y(i)-xold(j)^2)^2/(2*R)))*w(j);
        w(j)=exp(-((y(i)-xold(j)^2)^2/(2*0.1))); %R
    end
    
    %归一化
    w=w/sum(w);
    %w/sum（w）与k*w/sum（k*w）结果一模一样
    %(2*pi*R)^(-0.5)是常数
    %w(j),如果每次都重采样，每次w(j)都会被设为1/n,也是常数
    %所以可以将他们去掉
    
    %重采样
    %N<1/sum(wi^2),若不是每次都重采样的化，那么权重更新就要把w(j)乘上去
    %生成数组c
    c=zeros(1,n);
    c(1)=w(1);
    for j=2:n
        c(j)=c(j-1)+w(j);
    end
    
    %首先要重采样n个粒子，粒子数要跟之前相同
    for j = 1:n
        a = unifrnd(0,1);
        for k=1:n
            if(a<c(k))
                xnew(j)=xold(k);
                break;
            end
        end
    end
   
    xold = xnew;
    %权重都重新设为1/n, 只有重采样得时候有这一步
    
    for j=1:n
        w(j)=1/n;
    end
    %把每一步的后验概率期望赋值跟xplus
    xplus(i)=sum(xnew)/n;
    
end

plot(t,x,'r',t, xplus, 'b', 'LineWidth',2)

%y=x^2+R 似然概率是一个多峰分布, y=4, x可能为2，或者-2
%如果问题本身的性质就是强烈的非线性，比如多峰分布这种，粒子滤波也不行
%粒子滤波的计算速度有问题，而且占很多内存
    
    