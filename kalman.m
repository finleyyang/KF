%卡尔曼滤波
%X(K) = F*X(K-1)+Q
%Y(K) = H*X(K)+R

%生成一段时间t
t=0.1:0.01:1;
L=length(t);
%生成真实信号x，以及观测y
%首次初始化
x=zeros(1,L);
y=x;
for i=1:L
    x(i)=t(i)^2;
    y(i)=x(i)+normrnd(0,0.1);
    y(i)=x(i)+normrnd(0,0.1);
end
%生成信号完成

%滤波算法

%观测方程Y(K)=X(K)+R R~N(0,1)
%预测方程
%模型一，
%X(K)=X(K-1)+Q
%Y(K)=X(K)+R
%Q~N(0,1)
F1 = 1;
H1 = 1;
Q1 = 0.1;
R1 = 1;
%初始化X(k)+
Xplus1=zeros(1,L);

%设置一个初始值，假设Xplus1(1)~N(0.01, 0.01^2)
Xplus(1)=1;
Pplus1=0.01^2;
%%%%卡尔曼滤波
%X(K)minus=F*X(K-1)plus
%P(K)minus=F*P(K-1)plus*F'+Q
%K=P(K)minus*H'*inv(H*P(K)minus*H'+R)
%X(K)plus=X(K)minus+K(y(k)-H*X(K)minus)
%P(K)plus=(1-K*H)*P(K)minus

for i = 2:L
    %%%预测步
    Xminus1=F1*Xplus1(i-1);
    Pminus1=F1*Pplus1*F1'+Q1;
    %%%更新步
    K1=(Pminus1*H1')/(H1*Pminus1*H1'+R1);
    Xplus1(i)=Xminus1+K1*(y(i)-H1*Xminus1);
    Pplus1=(1-K1*H1)*Pminus1;
end

%plot(t,x,'r',t,y,'g',t,Xplus1,'b','LineWidth',2);
    
%模型二，
%X(K)=X(K-1)+X'(K-1)*dt+X''(K-1)*dt^2*(1/2!)+Q2
%Y(K)=X(K)+R R~N(0,1)
%Q~N(0,1)
%此时状态变量X=[X(K) X'(K) X''(K)]T(列向量)
%F = 1  dt  0.5dt^2
%    0, 1   dt
%    0, 0,  1
%H=[1, 0, 0]
%Q = Q2, 0, 0
%    0, Q3, 0
%    0, 0, Q4
dt=t(2)-t(1);
F2 = [1,dt,0.5*dt^2;0,1,dt;0,0,1];
H2 = [1, 0, 0];
Q2 = [0.2, 0, 0;0, 0.1, 0;0, 0, 0.01];
R2 = 5;

Xplus2=zeros(3, L);
Xplus2(1,1)=0.1^2;
Xplus2(2,1)=0;
Xplus2(3,1)=0;
Pplus2=[0.01, 0, 0;0, 0.01, 0; 0, 0, 0.001];

for i=2:L
    Xminus2=F2*Xplus2(:,i-1);
    Pminus2=F2*Pplus2*F2'+Q2;
    
    K2=(Pminus2*H2')/(H2*Pminus2*H2'+R2);
    Xplus2(:,i)=Xminus2+K2*(y(i)-H2*Xminus2);
    Pplus2=(eye(3)-K2*H2)*Pminus2;
end

plot(t,x,'r',t,y,'g',t,Xplus2(1,:),'b','LineWidth',2);



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



 