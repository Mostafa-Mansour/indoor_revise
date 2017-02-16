
%%%% MOdeling using Peter Croke toolbox %%%% 
close all
clear all

%% Robot Modeling 
rd=0.05;            %odometry RMS error in distance
rh=0.3*pi/180;      %odometry RMS error in heading
W=diag([rd,rh].^2); %Odometry error covariance matrix
veh=Vehicle(W);     % an object from vehicle class

%% Map Modeling
slength=10;
nlmarks=4;             %Number of the landmarks in the map
map=Map(nlmarks,slength);
figure
hold on
map.plot();
axis equal
%% Showing the initial robot pose in the map
pose=veh.x();
plot(pose(1),pose(2),'*');

%% Sensor Modeling
rr=0.1;             %RMS range sensor error
rb=pi/180;          % RMS bearing sensor error
V=diag([rr,rb].^2); % Sensor error Covariance Matrix
sensor=RangeBearingSensor(veh,map,V);

%% Measurment accusition
z=zeros(nlmarks,2);
for n=1:nlmarks
z(n,:)=sensor.h(veh.x()',n);
end
