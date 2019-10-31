% Kalman�˲��������ܽ��з���
%����������Kalman�˲����ܵ�Ӱ��
%% ������ʼ��
clc;
clear;
T = 1;%�״�ɨ������
N = 80/T;%�ܵĲ�������
M = 300;%Monte Carloģ�����
W = 20;%�ı�qs/rsģ�����
X = zeros(4,N,M);%Ŀ����ʵλ�á��ٶȡ�ģ�����
Z = zeros(2,N,M);%��������λ�õĹ۲�[x,y]
Result_Kal = zeros(1,W);%�ı�qs/rs��M��ʵ���ƽ��ͳ��

for k = 1:W
    qs = 0.01*k;
    Q1 = 0.01*[T.^3./3,T.^2./2,0,0;T.^2./2,T,0,0; 0,0,T.^3./3,T.^2./2;0,0,T.^2./2,T];%����������ֵ
    Q = qs*[T.^3./3,T.^2./2,0,0;T.^2./2,T,0,0; 0,0,T.^3./3,T.^2./2;0,0,T.^2./2,T];%����������ֵ
    %Q = qs*diag([1,1,1,1]) ;%����������ֵ
    rs = 10;
    R = rs*eye(2);%�۲�������ֵ
    F = [1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%״̬ת�ƾ���
    H = [1,0,0,0;0,0,1,0];%�۲����

    for i = 1:M
        X(:,1,i) = [0,1,0,1];%Ŀ���ʼλ�á��ٶ�[x,vx,y,vy]
        Z(:,1,i) = [X(1,1),X(3,1)];%�۲��ʼ��
        for t = 2:N
            %X(:,t,i) = F*X(:,t-1,i) + sqrtm(Q)*randn(4,1);
            %Z(:,t,i) = H*X(:,t,i) + sqrtm(R)*randn(2,1); 
            X(:,t,i) = F*X(:,t-1,i) + mvnrnd([0,0,0,0],Q1,1)';      
            Z(:,t,i) = H*X(:,t,i) + mvnrnd([0,0],R,1)';
        end
    end


    %% �������˲�
    Xkf = zeros(4,N,M);%Ŀ����ʵλ�á��ٶȡ�ģ�����
    Xkf(:,1,:) = X(:,1,:)+2;%�������˲�״̬��ʼ��
    P0 = eye(4);%Э�������ʼ��
    for j = 1:M
        for i = 2:N
            Xn = F*Xkf(:,i-1,j);%Ԥ��
            P1 = F*P0*F' + Q;%Ԥ�����Э����
            K = P1*H'/(H*P1*H' + R);%����
            Xkf(:,i,j) = Xn + K*(Z(:,i,j)-H*Xn);%״̬����
            P0 = (eye(4)-K*H)*P1;%�˲����Э�������
        end
    end
    %% �����¡��ۼ�
    Err_Pre = zeros(M,N);%ĳһ��ʵ���˲�ǰ��N�����ݵ����,M��ʵ��
    Err_Kal = zeros(M,N);%ĳһ��ʵ���˲���N�����ݵ����,M��ʵ��
    RMSE_Pre = zeros(2,N);%�˲�ǰM��ʵ��ľ�������x,y��ͬ����
    RMSE_Kal = zeros(2,N);%�˲���M��ʵ��ľ�������x,y��ͬ����
    for j = 1:N
        for i = 1:M
            Err_Pre(i,j) = EUD(X(:,j,i),Z(:,j,i));%�˲�ǰ�����ݵ������֮���ֱ�߾������
            Err_Kal(i,j) = EUD(X(:,j,i),Xkf(:,j,i));%�˲�������ݵ������֮���ֱ�߾������
        end
        RMSE_Pre(1,j) = sqrt( sum ((X(1,j,:)-Z(1,j,:)).^2)/M);%x����
        RMSE_Kal(1,j) = sqrt( sum ((X(1,j,:)-Xkf(1,j,:)).^2)/M);%x����
        RMSE_Pre(2,j) = sqrt (sum ((X(3,j,:)-Z(2,j,:)).^2)/M);%y����
        RMSE_Kal(2,j) = sqrt( sum ((X(3,j,:)-Xkf(3,j,:)).^2)/M);%y����             
    end
    Result_Kal(1,k) = sum(RMSE_Kal(1,20:N))/(N-20);%RMSE��̬ƽ��ֵ
    Result_Kal(2,k) = sum(RMSE_Kal(2,20:N))/(N-20);%RMSE��̬ƽ��ֵ
%     %% �������
%     [res_Pre,x_Pre] = hist(RMSE_Pre,M/2);%��Sum_RMSE_Pre�е���С�����ֵ֮�������ȷ�M/10�ݣ�
%     [res_Kal,x_Kal] = hist(RMSE_Kal,M/2);%���ص�result_Pre�����ڸ������ڵĸ�����x_Pre�Ǹ������������λ�����ꡣ
%     temp_Pre = 0; %�м�������˲�ǰ���Ƶ�Ⱥ�
%     temp_Kal = 0; %�м������Ƶ�Ⱥ����Ƶ�Ⱥ�
%     for i = 1:M/2
%         temp_Pre = temp_Pre + res_Pre(i);
%         temp_Kal = temp_Kal + res_Kal(i);
%         P_Pre(i) = temp_Pre.*100./(sum(res_Pre)+1);%�˲�ǰ���Ƶ�Ⱥ�-->Ƶ��
%         P_Kal(i) = temp_Kal.*100./(sum(res_Kal)+1);%�˲������Ƶ�Ⱥ�-->Ƶ��
%     end
 
    %% ��ͼ
    %ĳһ����켣ͼ
    if k == 1
        figure(1)
        hold on;box on;
        plot(X(1,:,1),X(3,:,1),'r-');  %��ʵ�켣
        plot(Z(1,:,1),Z(2,:,1),'bo');  %�����
        plot(Xkf(1,:,1),Xkf(3,:,1),'-ks','MarkerFace','g');  %�������˲��켣,%MarkerFace:���ݵ�ʵ�����
        xlabel('X����/m');
        ylabel('Y����/m');
        legend('��ʵ�켣','�����','�˲��켣')

        %ĳһ�������˲�N�����ݵ�ֱ�߾������ͼ
        figure(2)
        hold on; box on;
        plot(Err_Pre(1,:),'g');
        plot(Err_Kal(1,:),'r');
        xlabel('����ʱ��(t/s)');
        ylabel('ĳһ��ʵ��ŷʽ�������/m');
        legend('�˲�ǰ���','�˲������');
        
        %M��ʵ����������ͼ,x����
        figure(3)
        hold on; box on;
        plot(RMSE_Pre(1,:),'g');
        plot(RMSE_Kal(1,:),'r');        
        xlabel('����ʱ��(t/s)');
        ylabel('x����λ�þ��������/m');
        legend('�˲�ǰ���RMSE��300��ʵ�飩','�˲������RMSE��300��ʵ�飩');
        
        %M��ʵ����������ͼ,y����
        figure(4)
        hold on; box on;
        plot(RMSE_Pre(2,:),'g');
        plot(RMSE_Kal(2,:),'r');
        xlabel('����ʱ��(t/s)');
        ylabel('y����λ�þ��������/m');
        legend('�˲�ǰ���RMSE��300��ʵ�飩','�˲������RMSE��300��ʵ�飩');
        hold off;

%         %M��ʵ���˲������Ƶ�ʷֲ�ֱ��ͼ
%         figure
%         bar(x_Kal,res_Kal);
 
%         %M��ʵ���˲������Ƶ���ۼƷֲ�ͼ
%         figure
%         plot(x_Kal,P_Kal,'b');    
    end
        %M��ʵ����������ͼ,x����
        figure(5)
        plot(RMSE_Kal(1,:),'r');
        xlabel('����ʱ��(t/s)');
        ylabel('�˲���x����λ�þ��������/m');
        legend('��������ϵ��qsȡ[0.01-0.2]��RMSE');
        hold on; box on;
end

%W���ı��������/��������ϵ����RMSE����
figure(6)
plot(0.01*(2:W),Result_Kal(1,2:W));
xlabel('qs');
ylabel('�˲���x����RMSE��̬ƽ��ֵ/m');

figure(7)
plot(0.01*(2:W),Result_Kal(2,2:W));
xlabel('qs');
ylabel('�˲���y����RMSE��̬ƽ��ֵ/m');


%% ����ŷʽ�����Ӻ���
function dist = EUD(X1,X2);
    if length(X2)<=2
        dist = sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 );  %�˲�ǰ�����⣬X,Z
    else
        dist = sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 );  %�˲��������⣬X,Xkf
    end
end



