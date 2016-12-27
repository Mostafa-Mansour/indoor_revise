
rand('state',10);randn('state',10);
m=4;
realizations=50;
X_sys=zeros(3,m,realizations);
P_sys=zeros(3,3,m,realizations);
% X=zeros(2,m,realizations);
% P=zeros(2,2,m,realizations);
samples=100000;
a=0.25;
width_x=5;
length_y=5;
x_real=repmat(a,1,realizations)+repmat((width_x-2*a),1,realizations).*rand(1,realizations);
y_real=repmat(a,1,realizations)+repmat((length_y-2*a),1,realizations).*rand(1,realizations);
c_real=zeros(1,realizations)+repmat(0.5,1,realizations).*rand(1,realizations);
sigma2=zeros(3,m,realizations);
parfor itr=1:realizations
    disp(['realization number ' num2str(itr) ]);
    [X_sys(:,:,itr),P_sys(:,:,:,itr)]=getSystematic(x_real(itr),y_real(itr),c_real(itr),samples,width_x,length_y);
%     [X(:,:,itr),P(:,:,:,itr)]=getBootStrapIndoorKnownCourse(x_real(itr),y_real(itr),samples,width_x,length_y);
    sigma2(:,:,itr)=(X_sys(:,:,itr)-repmat([x_real(itr);y_real(itr);c_real(itr)],1,m));
    close all;
end

sigma=sqrt(sum(sigma2.^2,3)./(realizations-1));
P__=sum(P_sys,4)./(realizations);
% figure
% subplot(311)
% plot(linspace(1,m,m),sigma(1,:),'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
% subplot(312)
% plot(linspace(1,m,m),sigma(2,:),'r',linspace(1,m,m),sqrt(squeeze(P__(2,2,:))),'b');
% subplot(313)
% plot(linspace(1,m,m),sigma(3,:),'r',linspace(1,m,m),sqrt(squeeze(P__(3,3,:))),'b');


figure (1) 
hold on
subplot(321)
hold on
plot(linspace(1,m,m),abs(sigma(1,:)'-sqrt(squeeze(P__(1,1,:)))),'g');%,'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(323)
hold on
plot(linspace(1,m,m),abs(sigma(2,:)'-sqrt(squeeze(P__(2,2,:)))),'r');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(325)
hold on
plot(linspace(1,m,m),abs(sigma(3,:)'-sqrt(squeeze(P__(3,3,:)))),'r');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(322)
hold on
plot(linspace(1,m,m),(sigma(1,:)'-sqrt(squeeze(P__(1,1,:)))),'r');%,'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(324)
hold on
plot(linspace(1,m,m),(sigma(2,:)'-sqrt(squeeze(P__(2,2,:)))),'r');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(326)
hold on
plot(linspace(1,m,m),(sigma(3,:)'-sqrt(squeeze(P__(3,3,:)))),'r');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')