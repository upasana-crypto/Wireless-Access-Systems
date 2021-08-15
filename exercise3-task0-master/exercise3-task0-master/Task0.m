clear all;
close all;
clc;

%% given parameters
 
% room corners
roomCorners = [ 0 0 16 16 ; -1 3.4 3.4 -1 ];

% anchor positions
ancPos = [ 2.0325 4.0325 6.0325 8.0325 10.0325 12.0325 ; 0 2.4 0 2.4 0 2.4 ];

% starting position
p_start = [0; 1.2];

% final position
p_end = [14;1.2];
 
% number of samples between start and end
Nsamp = 189;

% discrete timestep duration between samples 
Dt = 70 * 10^(-3);

%% processing of the given data

% calculation of constant speed
 
distCov = sqrt(abs(p_start(1)-p_end(1))^2+abs(p_start(2)-p_end(2))^2);
timeInt = (Nsamp)*Dt;

speed = distCov/timeInt;
 
%% question b: 

% create time Axis
timeAx = zeros(1,Nsamp);
for i = 2:Nsamp
   timeAx(i) = i*Dt;
end

% first create the trajectory 
trajectory = zeros(2,Nsamp);
trajectory(2,:) = 1.2*ones(1,Nsamp);
for i = 1:Nsamp
    trajectory(1,i) = p_start(1,1)+timeAx(i)*speed;
end

% for each pos of the trajectory --> find d to all anchors
distanceToAnchors = zeros(Nsamp,size(ancPos,2));
for i = 1:Nsamp
    for j=1:size(ancPos,2)
        distanceToAnchors(i,j) = sqrt( abs( trajectory(1,i) - ancPos(1,j ) )^2  + abs( trajectory(2,i) - ancPos(2,j))^2 );
    end
end

%% question c: simulate RSSI values according to PLM

% path loss exponent a
a = 1.8;

% fixed distance d0
d0 = 0.3;

% reference RSSI at d0
RSSI_d0 = -38;

RSSI_SimNoNoise = zeros(size(distanceToAnchors));

for i =1:Nsamp
    for j=1:size(ancPos,2)
        RSSI_SimNoNoise(i,j)= RSSI_d0 - 10*a*log10(distanceToAnchors(i,j)/d0);
    end
end

%% question d: adding noise to simulated RSSI values 

% standart deviatiation
standDeviation_3dB =3;

% convert to decimal values
standDeviation_dec_3dB = 10^(standDeviation_3dB/10);

% generate Gaussian noise with given standart deviation
noise_3dB = standDeviation_dec_3dB.*randn(size(RSSI_SimNoNoise));

% add noise to the simulated RSSI values
RSSI_SimNoisy_3dB = RSSI_SimNoNoise+noise_3dB;

%% question e and f: write cost function and obtain pos estimates

% 1. obtain ML distance estimates
distanceHat_3dB = zeros(size(distanceToAnchors));
for i =1:Nsamp
    for j=1:size(ancPos,2)
        distanceHat_3dB(i,j)= d0*10^((-RSSI_SimNoisy_3dB(i,j)+RSSI_d0)/(10*a));
    end
end

posEstim_3dB = zeros(size(trajectory));
for i=1:Nsamp
    % iterate across all different samples
    
    % isolate the distances from ML for the current trajectory point
    distanceHat_vector_3dB = distanceHat_3dB(i,:);
    
    % 2. define cost function
    costFunction = @(posEst,distanceHat_vector,ancPos) log(distanceHat_vector)-log(sqrt(abs(posEst(1)-ancPos(1,:)).^2+abs(posEst(2)-ancPos(2,:)).^2));
   
    % 3. obtain position estimates
    posEstim_3dB(:,i) = lsqnonlin(@(posEst) costFunction(posEst,distanceHat_vector_3dB,ancPos),p_start);
end

% for evaluation purposes: calculate RMSE
calculate_RMSE=@(a,b) sqrt(mean((a(1,:)-b(1,:)).^2+(a(2,:)-b(2,:)).^2));
RMSE_3dB = calculate_RMSE(posEstim_3dB,trajectory);

%% question h : very noisy data

% standart deviatiation
standDeviation_10dB =10;

% convert to decimal values
standDeviation_dec_10dB = 10^(standDeviation_10dB/10);

% generate Gaussian noise with given standart deviation
noise_10dB = standDeviation_dec_10dB.*randn(size(RSSI_SimNoNoise));

% add noise to the simulated RSSI values
RSSI_SimNoisy_10dB = RSSI_SimNoNoise+noise_10dB;

% 1. obtain ML distance estimates
distanceHat_10dB = zeros(size(distanceToAnchors));
for i =1:Nsamp
    for j=1:size(ancPos,2)
        distanceHat_10dB(i,j)= d0*10^((-RSSI_SimNoisy_10dB(i,j)+RSSI_d0)/(10*a));
    end
end

posEstim_10dB = zeros(size(trajectory));
for i=1:Nsamp
    % iterate across all different samples
    
    % isolate the distances from ML for the current trajectory point
    distanceHat_vector_10dB = distanceHat_10dB(i,:);
    
    % 3. obtain position estimates
    posEstim_10dB(:,i) = lsqnonlin(@(posEst) costFunction(posEst,distanceHat_vector_10dB,ancPos),p_start);
end
RMSE_10dB = calculate_RMSE(posEstim_10dB,trajectory);

%% figures for questions g and h

figure;
subplot(2,1,1);
% 1. plot the room and the anchor positions
p1 = plot(roomCorners(1,:),roomCorners(2,:),'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% extend the plot limits to account for pos estimations that lie outside
x_axis_min = min(roomCorners(1,:))-2.5;
x_axis_max = max(roomCorners(1,:))+2.5;
y_axis_min = min(roomCorners(2,:))-2.5;
y_axis_max = max(roomCorners(2,:))+2.5;
xlim([x_axis_min x_axis_max]);
ylim([y_axis_min y_axis_max]);
% add a line to move the axis
line_x = roomCorners(1,1):roomCorners(1,end);
coord_y = min(roomCorners(2,:))*ones(size(line_x));
hold on;
plot(line_x,coord_y,'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the anchor positions
p2 = plot(ancPos(1,:),ancPos(2,:),'ro','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the true trajectory
p3 = plot(trajectory(1,:),trajectory(2,:),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the estimated trajectory for 3dB
p4 = plot(posEstim_3dB(1,:),posEstim_3dB(2,:),'g+','LineWidth', 3, 'MarkerSize',10);
hold on;
txt1 = ['RMSE-3dB = ' num2str(RMSE_3dB) ' m'];
my_leg_title = {txt1};
lgd = legend([p1 p2 p3 p4],'Room Frame','Anchor Positions','True trajectory','Est Pos 3dB');
lgd.FontSize = 12;
t = annotation('textbox',[0.7 0.55 0.1 0.1],'String',my_leg_title); 
sz = t.FontSize;
t.FontSize = 16;
title('Position Estimation with simulation and st-dev = 3dB','FontSize',16);
xlabel('x position');
ylabel('y position');

subplot(2,1,2);
% 1. plot the room and the anchor positions
p1 = plot(roomCorners(1,:),roomCorners(2,:),'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% extend the plot limits to account for pos estimations that lie outside
x_axis_min = min(roomCorners(1,:))-2.5;
x_axis_max = max(roomCorners(1,:))+2.5;
y_axis_min = min(roomCorners(2,:))-2.5;
y_axis_max = max(roomCorners(2,:))+2.5;
xlim([x_axis_min x_axis_max]);
ylim([y_axis_min y_axis_max]);
% add a line to move the axis
line_x = roomCorners(1,1):roomCorners(1,end);
coord_y = min(roomCorners(2,:))*ones(size(line_x));
hold on;
plot(line_x,coord_y,'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the anchor positions
p2 = plot(ancPos(1,:),ancPos(2,:),'ro','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the true trajectory
p3 = plot(trajectory(1,:),trajectory(2,:),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the estimated trajectory for 3dB
p4 = plot(posEstim_10dB(1,:),posEstim_10dB(2,:),'g+','LineWidth', 3, 'MarkerSize',10);
txt2 = ['RMSE-10dB = ' num2str(RMSE_10dB) ' m'];
my_leg_title = {txt2};
hold off;
lgd = legend([p1 p2 p3 p4],'Room Frame','Anchor Positions','True trajectory','Est Pos 10dB');
lgd.FontSize = 12;
t = annotation('textbox',[0.7 0.05 0.1 0.1],'String',my_leg_title); 
sz = t.FontSize;
t.FontSize = 16;
title('Position Estimation with simulation and st-dev = 10dB','FontSize',16);
xlabel('x position');
ylabel('y position');
