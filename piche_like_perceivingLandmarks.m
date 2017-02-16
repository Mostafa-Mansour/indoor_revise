function piche_like_perceivingLandmarks

close all
clear all
% Simulation and estimation of a nonlinear robot pose estimation problem
% There are some differences from the textbook model (Examples 5.3,7.1,7.2):

%% Define simulation parameters
lx=5;               % room width
ly=5;               % room length
r=0.25;             % measurement variance
x0=[lx/2;ly/2;0];   % initial state (x,y coordinates and heading)
DT=1;               % time step size
nk=500;             % number of steps
%mskip= 0;           % number of time steps per measurement

%% Define Landmarks

lmarks=[lx   lx     0       0;
           ly   0      0       ly];
       
X=zeros(2,nk);
figure(1)
% hold on 
grid on;
plot(lmarks(1,:),lmarks(2,:),'*','linewidth',20);

%% Measurement modeling 
x0=2;y0=1;
getRange=@(x,y) sqrt((lmarks(1,:)-x).^2 +(lmarks(2,:)-y).^2);

for t=1:nk
x=x0;y=y0;
x=x+randn();    y=y+randn();
realRange=getRange(x,y);
X(:,t)=[x;y];
plot(lmarks(1,:),lmarks(2,:),'*','linewidth',20);
hold on
plot([x lmarks(1,1)],[y lmarks(2,1)],'k--',[x lmarks(1,2)],[y lmarks(2,2)],'k--',...
     [x lmarks(1,3)],[y lmarks(2,3)],'k--',[x lmarks(1,4)],[y lmarks(2,4)],'k--');
 drawnow();
hold off
%  fprintf('Please press Enter \n');
 pause(0.5);
end


noise= makedist('Normal','sigma',sqrt(r));
plot
simTitle='robot pose estimation';

%% Simulations 
for isim=1:length(noise)
    disp(['Simulation ',num2str(isim)])
    rng('default'); rng(0)      % random number generator's starting value
    
    %% Track
    

% set seeding
% randn('state',0); rand('state',0);

%Motion model
model_f=@(x)x; %stationary mobile robot

%sesnor model
model_h=@(x,landmark) sqrt((landmark(1)-x(1))^2+(landmark(2)-x(2))^2);

%% Prior uniform distribution %%
a=0.25;
%%%%%%%%samples  %%%%%%%%%%%%%%%
samplePoints=repmat([a;a],1,samples)+repmat([(width_x-2*a);(length_y-2*a)],1,samples).*rand(2,samples);
xx=samplePoints(1,:);yy=samplePoints(2,:);

%% measurment modeling %%
% 8 measurements, 4 to every corner and 4 to the centers of the walls

landmarks=[width_x   width_x     0       0       ;%        width_x       width_x/2   0             width_x/2      width_x/4       3*width_x/4 ;
          length_y   0           0       length_y];%        length_y/2    0           length_y/2    length_y       length_y*3/4    length_y/4 ];

y=[];
rms_noise=0.05;
for n=1:length(landmarks)
    y(n)=model_h([x_real,y_real],landmarks(:,n))+rms_noise*randn();

end

%% Arrays for the results
X=zeros(2,length(landmarks));
P=zeros(2,2,length(landmarks));

%% Prior distibution 
prior=ones(1,samples);
prior=prior./samples;

%%Array for posterior distribution
posterior=zeros(length(landmarks)+1,samples);
posterior(1,:)=prior;
maxliklehood_values=zeros(1,samples);
cov__=zeros(2,2,samples);

%% loop for particle filter
figure
for n=1:length(landmarks) % loop for every meausrement (every land mark)
   for i=1:samples % loop for every particle
      sxi=model_h([xx(i),yy(i)],landmarks(:,n));
      maxliklehood_values(i)=Maximum_Likelihood_calculation_for_LRF(y(n),sxi,rms_noise^2);
      cov__(:,:,i)=[xx(i);yy(i)]*[xx(i);yy(i)]';
       
   end
   
   %% normalization
   posterior(n+1,:)=(maxliklehood_values.*posterior(n,:))./sum((maxliklehood_values.*posterior(n,:)));
    
   %% State vector Estimation %%
   X(1,n)=sum(xx.*posterior(n+1,:));
   X(2,n)=sum(yy.*posterior(n+1,:));
   
   %% Covariance Matrix
   for kk=1:samples
          cov__(:,:,kk)=posterior(n+1,kk).*cov__(:,:,kk);
                                    
   end               
   P(:,:,n)=sum(cov__,3)-[X(1,n);X(2,n)]*[X(1,n),X(2,n)];
   
  %% Resampling %%%
  [out]=ResamplingFast(samples,samples,posterior(n+1,:),[xx;yy]);
  xx=out(1,:);
  yy=out(2,:);
  
  %% Display results
  disp(['measurement  is ' num2str(n)]);
  plot(xx,yy,'.r');
  xlim([0 width_x])
  ylim([0 length_y])
  drawnow();
   
end
end
