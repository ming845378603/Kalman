% Kalman滤波估计性能进行分析
%过程噪声对Kalman滤波性能的影响
%% 参数初始化
clc;
clear;
T = 1;%雷达扫描周期
N = 80/T;%总的采样次数
M = 300;%Monte Carlo模拟次数
W = 20;%改变qs/rs模拟次数
X = zeros(4,N,M);%目标真实位置、速度、模拟次数
Z = zeros(2,N,M);%传感器对位置的观测[x,y]
Result_Kal = zeros(1,W);%改变qs/rs后，M次实验的平均统计

for k = 1:W
    qs = 0.01*k;
    Q1 = 0.01*[T.^3./3,T.^2./2,0,0;T.^2./2,T,0,0; 0,0,T.^3./3,T.^2./2;0,0,T.^2./2,T];%过程噪声均值
    Q = qs*[T.^3./3,T.^2./2,0,0;T.^2./2,T,0,0; 0,0,T.^3./3,T.^2./2;0,0,T.^2./2,T];%过程噪声均值
    %Q = qs*diag([1,1,1,1]) ;%过程噪声均值
    rs = 10;
    R = rs*eye(2);%观测噪声均值
    F = [1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%状态转移矩阵
    H = [1,0,0,0;0,0,1,0];%观测矩阵

    for i = 1:M
        X(:,1,i) = [0,1,0,1];%目标初始位置、速度[x,vx,y,vy]
        Z(:,1,i) = [X(1,1),X(3,1)];%观测初始化
        for t = 2:N
            %X(:,t,i) = F*X(:,t-1,i) + sqrtm(Q)*randn(4,1);
            %Z(:,t,i) = H*X(:,t,i) + sqrtm(R)*randn(2,1); 
            X(:,t,i) = F*X(:,t-1,i) + mvnrnd([0,0,0,0],Q1,1)';      
            Z(:,t,i) = H*X(:,t,i) + mvnrnd([0,0],R,1)';
        end
    end


    %% 卡尔曼滤波
    Xkf = zeros(4,N,M);%目标真实位置、速度、模拟次数
    Xkf(:,1,:) = X(:,1,:)+2;%卡尔曼滤波状态初始化
    P0 = eye(4);%协方差阵初始化
    for j = 1:M
        for i = 2:N
            Xn = F*Xkf(:,i-1,j);%预测
            P1 = F*P0*F' + Q;%预测误差协方差
            K = P1*H'/(H*P1*H' + R);%增益
            Xkf(:,i,j) = Xn + K*(Z(:,i,j)-H*Xn);%状态更新
            P0 = (eye(4)-K*H)*P1;%滤波误差协方差更新
        end
    end
    %% 误差更新、累计
    Err_Pre = zeros(M,N);%某一次实验滤波前，N个数据点误差,M次实验
    Err_Kal = zeros(M,N);%某一次实验滤波后，N个数据点误差,M次实验
    RMSE_Pre = zeros(2,N);%滤波前M次实验的均方根误差，x,y不同方向
    RMSE_Kal = zeros(2,N);%滤波后M次实验的均方根误差，x,y不同方向
    for j = 1:N
        for i = 1:M
            Err_Pre(i,j) = EUD(X(:,j,i),Z(:,j,i));%滤波前各数据点的误差，点之间的直线距离误差
            Err_Kal(i,j) = EUD(X(:,j,i),Xkf(:,j,i));%滤波后各数据点的误差，点之间的直线距离误差
        end
        RMSE_Pre(1,j) = sqrt( sum ((X(1,j,:)-Z(1,j,:)).^2)/M);%x方向
        RMSE_Kal(1,j) = sqrt( sum ((X(1,j,:)-Xkf(1,j,:)).^2)/M);%x方向
        RMSE_Pre(2,j) = sqrt (sum ((X(3,j,:)-Z(2,j,:)).^2)/M);%y方向
        RMSE_Kal(2,j) = sqrt( sum ((X(3,j,:)-Xkf(3,j,:)).^2)/M);%y方向             
    end
    Result_Kal(1,k) = sum(RMSE_Kal(1,20:N))/(N-20);%RMSE稳态平均值
    Result_Kal(2,k) = sum(RMSE_Kal(2,20:N))/(N-20);%RMSE稳态平均值
%     %% 概率求解
%     [res_Pre,x_Pre] = hist(RMSE_Pre,M/2);%将Sum_RMSE_Pre中的最小和最大值之间的区间等分M/10份，
%     [res_Kal,x_Kal] = hist(RMSE_Kal,M/2);%返回的result_Pre是落在该区间内的个数，x_Pre是该区间的中心线位置坐标。
%     temp_Pre = 0; %中间变量，滤波前误差频度和
%     temp_Kal = 0; %中间变量，频度和误差频度和
%     for i = 1:M/2
%         temp_Pre = temp_Pre + res_Pre(i);
%         temp_Kal = temp_Kal + res_Kal(i);
%         P_Pre(i) = temp_Pre.*100./(sum(res_Pre)+1);%滤波前误差频度和-->频率
%         P_Kal(i) = temp_Kal.*100./(sum(res_Kal)+1);%滤波后误差频度和-->频率
%     end
 
    %% 画图
    %某一试验轨迹图
    if k == 1
        figure(1)
        hold on;box on;
        plot(X(1,:,1),X(3,:,1),'r-');  %真实轨迹
        plot(Z(1,:,1),Z(2,:,1),'bo');  %量测点
        plot(Xkf(1,:,1),Xkf(3,:,1),'-ks','MarkerFace','g');  %卡尔曼滤波轨迹,%MarkerFace:数据点实心填充
        xlabel('X方向/m');
        ylabel('Y方向/m');
        legend('真实轨迹','量测点','滤波轨迹')

        %某一次试验滤波N个数据点直线距离误差图
        figure(2)
        hold on; box on;
        plot(Err_Pre(1,:),'g');
        plot(Err_Kal(1,:),'r');
        xlabel('跟踪时刻(t/s)');
        ylabel('某一次实验欧式距离误差/m');
        legend('滤波前误差','滤波后误差');
        
        %M次实验均方根误差图,x方向
        figure(3)
        hold on; box on;
        plot(RMSE_Pre(1,:),'g');
        plot(RMSE_Kal(1,:),'r');        
        xlabel('跟踪时刻(t/s)');
        ylabel('x方向位置均方根误差/m');
        legend('滤波前误差RMSE（300次实验）','滤波后误差RMSE（300次实验）');
        
        %M次实验均方根误差图,y方向
        figure(4)
        hold on; box on;
        plot(RMSE_Pre(2,:),'g');
        plot(RMSE_Kal(2,:),'r');
        xlabel('跟踪时刻(t/s)');
        ylabel('y方向位置均方根误差/m');
        legend('滤波前误差RMSE（300次实验）','滤波后误差RMSE（300次实验）');
        hold off;

%         %M次实验滤波后误差频率分布直方图
%         figure
%         bar(x_Kal,res_Kal);
 
%         %M次实验滤波后误差频率累计分布图
%         figure
%         plot(x_Kal,P_Kal,'b');    
    end
        %M次实验均方根误差图,x方向
        figure(5)
        plot(RMSE_Kal(1,:),'r');
        xlabel('跟踪时刻(t/s)');
        ylabel('滤波后x方向位置均方根误差/m');
        legend('过程噪声系数qs取[0.01-0.2]的RMSE');
        hold on; box on;
end

%W个改变过程噪声/量测噪声系数的RMSE分析
figure(6)
plot(0.01*(2:W),Result_Kal(1,2:W));
xlabel('qs');
ylabel('滤波后x方向RMSE稳态平均值/m');

figure(7)
plot(0.01*(2:W),Result_Kal(2,2:W));
xlabel('qs');
ylabel('滤波后y方向RMSE稳态平均值/m');


%% 计算欧式距离子函数
function dist = EUD(X1,X2);
    if length(X2)<=2
        dist = sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );  %滤波前误差求解，X,Z
    else
        dist = sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );  %滤波后误差求解，X,Xkf
    end
end



