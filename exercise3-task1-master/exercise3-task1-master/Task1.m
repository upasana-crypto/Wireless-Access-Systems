clear all;
close all;
clc;

%% load file
trajectory_straight = load('trajectory_straight.mat');

%% exctract data

% true position
posTrue = trajectory_straight.Data.posTrue;

% anchor positions
ancPos = trajectory_straight.Data.ancPos;

% measured RSSI values
RSSI_meas = trajectory_straight.Data.RSSI_meas;

% room coordinates
room = trajectory_straight.Data.room;

%% data processing

% define Number of samples
Nsamp = size(posTrue,2);

% calculate distance from anchors for each pos of the trajectory
distanceToAnchors = zeros(Nsamp,size(ancPos,2));
for i = 1:Nsamp
    for j=1:size(ancPos,2)
        distanceToAnchors(i,j) = sqrt( abs( posTrue(1,i) - ancPos(1,j ) )^2  + abs( posTrue(2,i) - ancPos(2,j))^2 );
    end
end

%% simulate RSSI values 

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

%% plot

figure;

% plot the measured RSSI values over distance for anchor 1
p1= plot(distanceToAnchors(:,1),RSSI_meas(:,1),'bo','LineWidth', 3, 'MarkerSize',10);
hold on;

% plot the simulated RSSI values over distance for anchor 1
p2 = plot(distanceToAnchors(:,1),RSSI_SimNoNoise(:,1),'k*','LineWidth', 3, 'MarkerSize',10);
hold on;

lgd = legend([p1 p2 ],'measured RSSI','simulated RSSI no-noise');
lgd.FontSize = 16;

title('Plot of RSSI values over distance for anchor 1','FontSize',16);
xlabel('distance from anchor 1 in meters');
ylabel('RSSI values in dB');

%% use code for position estimation from Task0

% 1. obtain ML distance estimates
distanceHat = zeros(size(distanceToAnchors));
for i =1:Nsamp
    for j=1:size(ancPos,2)
        distanceHat(i,j)= d0*10^((-RSSI_meas(i,j)+RSSI_d0)/(10*a));
    end
end

posEstim = zeros(size(posTrue));
for i=1:Nsamp
    % iterate across all different samples
    
    % isolate the distances from ML for the current trajectory point
    distanceHat_vector = distanceHat(i,:);
    
    % 2. define cost function
    costFunction = @(posEst,distanceHat_vector,ancPos) log(distanceHat_vector)-log(sqrt(abs(posEst(1)-ancPos(1,:)).^2+abs(posEst(2)-ancPos(2,:)).^2));

    % 3. obtain position estimates
    posEstim(:,i) = lsqnonlin(@(posEst) costFunction(posEst,distanceHat_vector,ancPos),posTrue(:,1));
end

%% for evaluation purposes: calculate RMSE
calculate_RMSE=@(a,b) sqrt(mean((a(1,:)-b(1,:)).^2+(a(2,:)-b(2,:)).^2));
RMSE = calculate_RMSE(posEstim,posTrue);

%% visualize localization before calibration
figure;
% 1. plot the room and the anchor positions
p1 = plot(room(1,:),room(2,:),'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% extend the plot limits to account for pos estimations that lie outside
x_axis_min = min(room(1,:))-2.5;
x_axis_max = max(room(1,:))+2.5;
y_axis_min = min(room(2,:))-2.5;
y_axis_max = max(room(2,:))+2.5;
xlim([x_axis_min x_axis_max]);
ylim([y_axis_min y_axis_max]);
% add a line to move the axis
line_x = room(1,1):room(1,end);
coord_y = min(room(2,:))*ones(size(line_x));
plot(line_x,coord_y,'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the anchor positions
p2 = plot(ancPos(1,:),ancPos(2,:),'ro','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the true trajectory
p3 = plot(posTrue(1,:),posTrue(2,:),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the estimated trajectory for 3dB
p4 = plot(posEstim(1,:),posEstim(2,:),'g+','LineWidth', 3, 'MarkerSize',10);
hold on;
txt1 = ['RMSE = ' num2str(RMSE) ' m'];
my_leg_title = {txt1};
hold off;
lgd = legend([p1 p2 p3 p4],'Room Frame','Anchor Positions','True position','Est Pos');
lgd.FontSize = 16;
t = annotation('textbox',[0.5 0.1 0.1 0.1],'String',my_leg_title); 
sz = t.FontSize;
t.FontSize = 16;
title('Position Estimation with measured RSSI no-calibration','FontSize',16);
xlabel('x position');
ylabel('y position');

%% define Ak matrix and perform the calibration
calibration_RSSId0_a = zeros(size(ancPos,2),2);
for j=1:size(ancPos,2)
    Ak = zeros(Nsamp,2);
    for i =1:Nsamp 
        Ak(i,1) = 1;
        Ak(i,2) = -10*log10( distanceToAnchors(i,j)/d0);
    end
    calibration_RSSId0_a(j,:) = Ak\RSSI_meas(:,j);
end

%% construct the updated PLM per anchor
RSSI_upd_PLM = zeros(size(RSSI_meas));
for j=1:size(ancPos,2)
    for i=1:Nsamp
        RSSI_upd_PLM(i,j) = calibration_RSSId0_a(j,1) - 10*log10(distanceToAnchors(i,j)/d0)*calibration_RSSId0_a(j,2);        
    end    
end

%% question g: calibration dataset as a scatterplot
figure
subplot(3,2,1)
% anchor 1
p1 = plot(distanceToAnchors(:,1),RSSI_meas(:,1),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
p2 = plot(distanceToAnchors(:,1),RSSI_upd_PLM(:,1),'k+','LineWidth', 3, 'MarkerSize',10);
title('scatterplot for RSSI values for anchor 1','FontSize',16);
lgnd = legend([p1 p2],'measured RSSI','PLM RSSI calibr','Location','NorthEast');
lgnd.FontSize = 16;

subplot(3,2,2)
% anchor 2
p1 = plot(distanceToAnchors(:,2),RSSI_meas(:,2),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
p2 = plot(distanceToAnchors(:,2),RSSI_upd_PLM(:,2),'k+','LineWidth', 3, 'MarkerSize',10);
title('scatterplot for RSSI values for anchor 2','FontSize',16);
lgnd = legend([p1 p2],'measured RSSI','PLM RSSI calibr','Location','NorthEast');
lgnd.FontSize = 16;

subplot(3,2,3)
% anchor 3
p1 = plot(distanceToAnchors(:,3),RSSI_meas(:,3),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
p2 = plot(distanceToAnchors(:,3),RSSI_upd_PLM(:,3),'k+','LineWidth', 3, 'MarkerSize',10);
title('scatterplot for RSSI values for anchor 3','FontSize',16);
lgnd = legend([p1 p2],'measured RSSI','PLM RSSI calibr','Location','NorthEast');
lgnd.FontSize = 16;

subplot(3,2,4)
% anchor 4
p1 = plot(distanceToAnchors(:,4),RSSI_meas(:,4),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
p2 = plot(distanceToAnchors(:,4),RSSI_upd_PLM(:,4),'k+','LineWidth', 3, 'MarkerSize',10);
title('scatterplot for RSSI values for anchor 4','FontSize',16);
lgnd = legend([p1 p2],'measured RSSI','PLM RSSI calibr','Location','NorthEast');
lgnd.FontSize = 16;

subplot(3,2,5)
% anchor 5
p1 = plot(distanceToAnchors(:,5),RSSI_meas(:,5),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
p2 = plot(distanceToAnchors(:,5),RSSI_upd_PLM(:,5),'k+','LineWidth', 3, 'MarkerSize',10);
title('scatterplot for RSSI values for anchor 5','FontSize',16);
lgnd = legend([p1 p2],'measured RSSI','PLM RSSI calibr','Location','NorthEast');
lgnd.FontSize = 16;

subplot(3,2,6)
% anchor 6
p1 = plot(distanceToAnchors(:,6),RSSI_meas(:,6),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
p2 = plot(distanceToAnchors(:,6),RSSI_upd_PLM(:,6),'k+','LineWidth', 3, 'MarkerSize',10);
title('scatterplot for RSSI values for anchor 5','FontSize',16);
lgnd = legend([p1 p2],'measured RSSI','PLM RSSI calibr','Location','NorthEast');
lgnd.FontSize = 16;

%% question j: estimate user positions from dataset trajectory_2.mat

% 1. load trajectory_2.mat
trajectory_2 = load('trajectory_2.mat');

% 2. extract RSSI values
RSSI_meas_Tr2 = trajectory_2.Data_2.RSSI_meas;

% 3. extract actual position
posTrue_Tr2 = trajectory_2.Data_2.posTrue;

% 3. obtain ML distance estimates using calibrated values
distHatTr2_calibr = zeros(size(RSSI_meas_Tr2));
for i =1:size(distHatTr2_calibr,1)
    for j=1:size(distHatTr2_calibr,2)
        distHatTr2_calibr(i,j)= d0*10^((-RSSI_meas_Tr2(i,j)+calibration_RSSId0_a(j,1))/(10*calibration_RSSId0_a(j,2)));
    end
end

% 4. use cost function and lsqnonlin to obtain the pos estimates
posEstim_Tr2_calibr = zeros(size(posTrue_Tr2));
for i=1:size(posEstim_Tr2_calibr,2)
    % iterate across all different samples
    
    % isolate the distances from ML for the current trajectory point
    distHatTr2_calibr_vec = distHatTr2_calibr(i,:);
    
    % obtain position estimates
    posEstim_Tr2_calibr(:,i) = lsqnonlin(@(posEst) costFunction(posEst,distHatTr2_calibr_vec,ancPos),posTrue_Tr2(:,1));
end

% 5. evaluate calibrated estimation of trajectory_2 using RMSE
RMSE_Tr2_calib = calculate_RMSE(posEstim_Tr2_calibr,posTrue_Tr2);

% 6. obtain ML distance estimates using Uncalibrated values
distHatTr2_uncalibr = zeros(size(RSSI_meas_Tr2));
for i =1:size(distHatTr2_uncalibr,1)
    for j=1:size(distHatTr2_uncalibr,2)
        distHatTr2_uncalibr(i,j)= d0*10^((-RSSI_meas_Tr2(i,j)+RSSI_d0)/(10*a));
    end
end

% 7. use cost function and lsqnonlin to obtain the pos estimates
posEstim_Tr2_uncalibr = zeros(size(posTrue_Tr2));
for i=1:size(posEstim_Tr2_uncalibr,2)
    % iterate across all different samples
    
    % isolate the distances from ML for the current trajectory point
    distHatTr2_uncalibr_vec = distHatTr2_uncalibr(i,:);
    
    % obtain position estimates
    posEstim_Tr2_uncalibr(:,i) = lsqnonlin(@(posEst) costFunction(posEst,distHatTr2_uncalibr_vec,ancPos),posTrue_Tr2(:,1));
end

% 8. evaluate uncalibrated estimation of trajectory_2 using RMSE
RMSE_Tr2_uncalib = calculate_RMSE(posEstim_Tr2_uncalibr,posTrue_Tr2);
%% visualize question j
figure;
subplot(2,1,1);
% 1. plot the room and the anchor positions
p1 = plot(room(1,:),room(2,:),'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% extend the plot limits to account for pos estimations that lie outside
x_axis_min = min(room(1,:))-2.5;
x_axis_max = max(room(1,:))+2.5;
y_axis_min = min(room(2,:))-2.5;
y_axis_max = max(room(2,:))+2.5;
xlim([x_axis_min x_axis_max]);
ylim([y_axis_min y_axis_max]);
% add a line to move the axis
line_x = room(1,1):room(1,end);
coord_y = min(room(2,:))*ones(size(line_x));
plot(line_x,coord_y,'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the anchor positions
p2 = plot(ancPos(1,:),ancPos(2,:),'ro','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the true trajectory
p3 = plot(posTrue_Tr2(1,:),posTrue_Tr2(2,:),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the estimated trajectory with calibration
p4 = plot(posEstim_Tr2_calibr(1,:),posEstim_Tr2_calibr(2,:),'g+','LineWidth', 3, 'MarkerSize',10);
hold on;
txt1 = ['RMSE = ' num2str(RMSE_Tr2_calib) ' m'];
my_leg_title = {txt1};
lgd = legend([p1 p2 p3 p4],'Room Frame','Anchor Positions','True position','Est Pos with calibr');
lgd.FontSize = 14;
t = annotation('textbox',[0.15 0.8 0.1 0.1],'String',my_leg_title); 
sz = t.FontSize;
t.FontSize = 16;
title('Position Estimation for trajectory 2 using calibration','FontSize',16);
xlabel('x position');
ylabel('y position');

subplot(2,1,2);
% 1. plot the room and the anchor positions
p1 = plot(room(1,:),room(2,:),'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% extend the plot limits to account for pos estimations that lie outside
x_axis_min = min(room(1,:))-2.5;
x_axis_max = max(room(1,:))+2.5;
y_axis_min = min(room(2,:))-2.5;
y_axis_max = max(room(2,:))+2.5;
xlim([x_axis_min x_axis_max]);
ylim([y_axis_min y_axis_max]);
% add a line to move the axis
line_x = room(1,1):room(1,end);
coord_y = min(room(2,:))*ones(size(line_x));
plot(line_x,coord_y,'k','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the anchor positions
p2 = plot(ancPos(1,:),ancPos(2,:),'ro','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the true trajectory
p3 = plot(posTrue_Tr2(1,:),posTrue_Tr2(2,:),'b*','LineWidth', 3, 'MarkerSize',10);
hold on;
% plot the estimated trajectory with calibration
p4 = plot(posEstim_Tr2_uncalibr(1,:),posEstim_Tr2_uncalibr(2,:),'g+','LineWidth', 3, 'MarkerSize',10);
hold on;
txt1 = ['RMSE = ' num2str(RMSE_Tr2_uncalib) ' m'];
my_leg_title = {txt1};
lgd = legend([p1 p2 p3 p4],'Room Frame','Anchor Positions','True position','Est Pos without calibr');
lgd.FontSize = 14;
t = annotation('textbox',[0.15 0.3 0.1 0.1],'String',my_leg_title); 
sz = t.FontSize;
t.FontSize = 16;
title('Position Estimation for trajectory 2 without calibration','FontSize',16);
xlabel('x position');
ylabel('y position');
