
m=4;
realizations=1000;
X=zeros(2,m,realizations);
P=zeros(2,2,m,realizations);
samples=100000;
a=0.25;
width_x=5;
length_y=5;
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
figure
subplot(211)
plot(linspace(1,m,m),sigma(1,:),'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
subplot(212)
plot(linspace(1,m,m),sigma(2,:),'r',linspace(1,m,m),sqrt(squeeze(P__(2,2,:))),'b');

figure (2)
hold on
subplot(231)
plot(linspace(1,m,m),abs(sigma(1,:)'-sqrt(squeeze(P__(1,1,:)))),'b');%,'r',linspace(1,m,m),sqrt(squeeze(P__(1,1,:))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
subplot(234)
plot(linspace(1,m,m),abs(sigma(2,:)'-sqrt(squeeze(P__(2,2,:)))),'b');
grid on,
xlabel('Number of Meaurements')
ylabel('Difference')
% subplot(313)
% plot(linspace(1,m,m),abs(sigma(3,:)'-sqrt(squeeze(P__(3,3,:)))),'b');