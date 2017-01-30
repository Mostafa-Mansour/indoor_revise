rand('state',0);randn('state',0);

m=4;
realizations=100;
X=zeros(2,m,realizations);
P=zeros(2,2,m,realizations);
samples=30000;
a=0.25;
width_x=3;
length_y=3;
x_real=repmat(a,1,realizations)+repmat((width_x-2*a),1,realizations).*rand(1,realizations);
y_real=repmat(a,1,realizations)+repmat((length_y-2*a),1,realizations).*rand(1,realizations);
sigma2=zeros(2,m,realizations);
parfor itr=1:realizations
    disp(['realization number ' num2str(itr) ]);
    [X(:,:,itr),P(:,:,:,itr)]=getBootStrapIndoorKnownCourse(x_real(itr),y_real(itr),samples,width_x,length_y);
    sigma2(:,:,itr)=(X(:,:,itr)-repmat([x_real(itr);y_real(itr)],1,m));
    close all;
end

sigma=sqrt(sum(sigma2.^2,3)./(realizations-1));
P__=sum(P,4)./realizations;
figure (7)
subplot(233)
hold on
plot(linspace(1,m,m),sigma(1,:),'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
subplot(236)
hold on
plot(linspace(1,m,m),sigma(2,:),'r',linspace(1,m,m),sqrt(squeeze(P__(2,2,:))),'b');

figure(4) 
hold on
subplot(232)
plot(linspace(1,m,m),abs(sigma(1,:)'-sqrt(squeeze(P__(1,1,:)))),'b');%,'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(235)
plot(linspace(1,m,m),abs(sigma(2,:)'-sqrt(squeeze(P__(2,2,:)))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
% subplot(313)
% plot(linspace(1,m,m),abs(sigma(3,:)'-sqrt(squeeze(P__(3,3,:)))),'b');