function [ X,P ] = getSystematic( x_real,y_real,c_real,samples,width_x,length_y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%%% Boot strap filter fot indoor navigation%%

 rand('state',0);randn('state',0);

%%%%%%%% The same parameters that used in kalman filter code %%%%%%%%%%%%
%%%%%%%%%Program constatnts %%%%%%%%%%%
% samples=100000;  % number of samples

a=0.25;
%%%%%%%%Room Dimensions %%%%%%%%%%%
% width_x=5;    % width of the room
% length_y=5;   % length of the room 

%%%%%%%%samples  %%%%%%%%%%%%%%%
samplePoints=repmat([a;a;0],1,samples)+repmat([(width_x-2*a);(length_y-2*a);0.5],1,samples).*rand(3,samples);
%%%%%%%%real values%%%%%%%%%%
% x_real=1;
% y_real=4;
% c_real=0.8;

%%%%%%%%Measurement modeling %%%%%%%
k=0;          % course as a constant   
           % number of measurements
           %[360-atand(x_real/(length_y-y_real)),
meas_angles=[atand((width_x-x_real)/(length_y-y_real)),360-atand(x_real/(length_y-y_real)),180-atand((width_x-x_real)/y_real),180+atand(x_real/y_real)];%,10,100,190,280];%,90,160,80,269,200,310];     % measurements angels
m=length(meas_angles);

real_range=[];
% real_range(:)=getRange(x_real,y_real,k+meas_angles(:),getWallNum(x_real,y_real,k+meas_angles(:),[width_x;length_y]),[width_x;length_y]);
for itr=1:m
    real_range(itr)=getRange_(x_real,y_real,k+meas_angles(itr),getWallNum(x_real,y_real,k+meas_angles(itr),[width_x;length_y]),[width_x;length_y]);
    
end




%%%%%Array for the results%%%%%
X=zeros(3,m);
% X(2,1,:)=y_real;
P=zeros(3,3,m);

%%%% noisy measurement %%%%%%% 
variance=.05^2;                   %variance of measurement white noise (5 cm )
R=variance*eye(3*m);
%%% prioir distribution %%%
prior=ones(1,samples);
prior=prior./samples;

%%%posterior distribution%%
posterior=zeros(m+1,samples);
posterior(1,:)=prior;
maxliklehood_values=zeros(1,samples);
cov__=zeros(3,3,samples);

%measurement modeling%
Y=real_range+c_real+sqrt(variance)*randn(1,length(real_range));

figure
for ii=1:m
    

for i=1:samples
                 
                sx=getRange_(samplePoints(1,i),samplePoints(2,i),k+meas_angles(ii),getWallNum(samplePoints(1,i),samplePoints(2,i),k+meas_angles(ii),[width_x;length_y]),[width_x;length_y]);
%                               
%                
                maxliklehood_values(i)=Maximum_Likelihood_calculation_for_LRF(Y(ii),sx+samplePoints(3,i),R(ii,ii));
      
                cov__(:,:,i)=samplePoints(:,i)*samplePoints(:,i)';
    
end
                sum_=sum(maxliklehood_values);
              if (sum_~=0)
                  
               %%%normalization%%%
               posterior(ii+1,:)=(maxliklehood_values.*posterior(ii,:))./sum((maxliklehood_values.*posterior(ii,:)));

               
               
               %%% State Vector Estimation %%%
               X(1,ii)=sum(samplePoints(1,:).*posterior(ii+1,:));
               X(2,ii)=sum(samplePoints(2,:).*posterior(ii+1,:));
               X(3,ii)=sum(samplePoints(3,:).*posterior(ii+1,:));
               
               %%%Covariance Matrix Estimation%%%
%                cov__=[squeeze(cov__(1,1,:)).*posterior(1,:)' squeeze(cov__(1,2,:)).*posterior(1,:);squeeze(cov__(2,1,:)).*posterior(1,:) squeeze(cov__(2,2,:)).*posterior(1,:)];
%                cov__=reshape(reshape(cov__,[1,40000]).*repmat(posterior(ii+1,:),1,4),[2,2,10000]);
               for kk=1:samples
                   cov__(:,:,kk)=posterior(ii+1,kk).*cov__(:,:,kk);
                                    
               end               
               P(:,:,ii)=sum(cov__,3)-[X(1,ii);X(2,ii);X(3,ii)]*[X(1,ii),X(2,ii),X(3,ii)];
                  
                %% Resampling %%%
               [out]=ResamplingFast(samples,samples,posterior(ii+1,:),samplePoints);
               samplePoints(1,:)=out(1,:);
               samplePoints(2,:)=out(2,:);
               samplePoints(3,:)=out(3,:);
              
%               else
%                   %%% State Vector Estimation %%%
%                X(1,ii)=sum(samplePoints(1,:).*posterior(ii+1,:));
%                X(2,ii)=sum(samplePoints(2,:).*posterior(ii+1,:));
%                
%                %%%Covariance Matrix Estimation%%%
%                for k=1:samples
%                    cov__(:,:,k)=posterior(ii+1,k).*cov__(:,:,k);
%                                     
%                end               
%                P(:,:,ii)=sum(cov__,3)-[X(1,ii);X(2,ii)]*[X(1,ii),X(2,ii)];
              end
              


disp(['measurement  is ' num2str(ii)]);
  plot(samplePoints(1,:),samplePoints(2,:),'.r');
  xlim([0 width_x])
  ylim([0 length_y])
  drawnow();
end

plot(squeeze(X(1,m,end)),squeeze(X(2,m,end)),'.r');



end

