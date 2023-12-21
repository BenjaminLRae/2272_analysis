% 
%  Laser tracker import and processing script
%  Author: Ben Rae
%  Last edited: 13/12/2023
%

%clc
clear
close all

% Import data using automated importer
Data = importfile('P:\NIN Projects\NIN2272\Laser Tracker Data\NIN2272\TEST_1_081223\TEST_1_M.txt');

% Number of points captured in file
numPoints = height(Data);

% Capture rate (Hz)
captureRate = 1000;

% Time between each point (s)
dt = 1 / captureRate;

% Total capture time (s)
captureTime = numPoints * dt;
captureTimeMin = captureTime / 60;

%% Calculate time and convert from inches for each point
disp('Calculating, please wait...');
for i=1:height(Data)
    Data.time(i) = i * dt; 
end

Data.x = Data.x * 25.4;
Data.y = Data.y * 25.4;
Data.z = Data.z * 25.4;

disp('Calculation complete');

%% Create normalised option (no negative axis moves)
DataNorm = Data;
DataNorm.x = abs(DataNorm.x);
DataNorm.y = abs(DataNorm.y);
DataNorm.z = abs(DataNorm.z);

% DataNorm.x = DataNorm.x;
% DataNorm.y = DataNorm.y;
% DataNorm.z = DataNorm.z;

%% Graphing
t = tiledlayout(3,6);
ax1 = nexttile(1, [3 3]);
hold on;
f3d = plot3(DataNorm.x, DataNorm.y, DataNorm.z);
xlim padded;
ylim padded;
zlim padded;
f3dOrigin = plot3(0,0,0);
f3dOrigin.Color = "red";
f3dOrigin.Marker = "o";
daspect([1 1 1]);
grid on;
title('3D view');
ylabel('distance (')

ax2 = nexttile(4, [1 3]);
hold on;
plot(DataNorm.time, DataNorm.x);
title('X axis');
grid on;

ax3 = nexttile(10, [1 3]);
hold on;
plot(DataNorm.time, DataNorm.y, 'red');
title('Y axis');
grid on;

ax4 = nexttile(16, [1 3]);
hold on;
plot(DataNorm.time, DataNorm.z, 'green');
title('Z axis');
grid on;
xlabel("time (s)");

%% peak detection
figPeak = figure;
title("Peak detection on X axis");

plot(DataNorm.time, DataNorm.x);

speed = gradient(DataNorm.x); % 1st derivative / speed
acc = gradient(speed); % 2nd derivative / acceleration 
smoothedAcc = smoothdata(acc,"movmean",400); % smoothing 2nd deriv. data to see peaks
smoothedAccAbs = abs(smoothedAcc); % make all points positive so that peak finding sees negative jerk too

figure 
plot(DataNorm.time, speed, 'b')
hold on
plot(DataNorm.time, acc, 'r')
plot(DataNorm.time, smoothedAcc, 'g')
hold off

% [pks,locs,p] = findpeaks(smoothedAccAbs, DataNorm.time,'MinPeakProminence',0.0001,'MinPeakHeight',1e-05);
% for i=1:length(locs)
%     xline(locs(i),'g');
% end

%% Speed Calculation

MyData = table2array(DataNorm);

Mydist=zeros(length(MyData),1);
speed=zeros(length(MyData),1);

% Taking the 1st derivative of position wrt time
time_array=0:dt:dt*(length(MyData)-1);
counter=1;
for p=1:length(MyData)-1
    P1=MyData(p,1:3);
    P2=MyData(p+1,1:3);
    Mydist(counter)=norm(P2-P1);
    speed(counter)=Mydist(counter)/dt*60; %convert inch to mm
    counter=counter+1;
end

% Smoothing the feed rate data
SmoothSpeed = smoothdata(speed);

% Feed rate plot
figure
plot(time_array,SmoothSpeed,'r','LineWidth',2)
title('Attained Feed Rate (mm/min)')
xlabel('Time (sec)')
ylabel('Feed Rate (mm/min)')
hold on 
plot(time_array,speed,'b')
% ylim([0 5000])
hold off

%% Acceleration calculation

Mydist=zeros(length(MyData),1);
acc=zeros(length(MyData),1);

% Taking the 2nd derivative of position wrt time - Smoothed feedrate
time_array=0:dt:dt*(length(speed)-1);
ind=1;
for p=1:length(speed)-1
    V1=speed(p);
    V2=speed(p+1);
    Mydist(ind)=V2-V1;
    acc(ind)=Mydist(ind)/dt/1000;
    ind=ind+1;
end

figure
plot(time_array,acc,'r')
title('Attained Acceleration (m/s^2)')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')

%%
jerk = gradient(acc);
smoothedJerk = smoothdata(jerk, "gaussian", 1000);
figure;
plot(DataNorm.time, smoothedJerk);
yyaxis right
plot(DataNorm.time, DataNorm.x);
title("Smoothed Jerk over Position");

figure;
plot(DataNorm.time, smoothedJerk);
yyaxis right
plot(DataNorm.time, smoothedAcc);
title("Smoothed Jerk over Acc");

%% Try to separate the movements based on acceleration. When acceleration is not zero, mark cut in data.

threshold = 0.0001; % difference from zero that implies acceleration/movement
inMove = false; % flag indicating whether we're in a move, i.e. not static

for i=0:length(grad) 
    if grad(i) > threshold
        % postive threshold exceeded
        if inMove
            % we're already in a move, do nothing.

        else 
            % we're not in a move, so set the inMove flag, and record the
            % index as the start of a move.
            inMove = true;
        end        
    elseif grad(i) < (0-threshold)
        % negative threshold exceeded
        if inMove
            % we're already in a move, do nothing.
        else 
            % we're not in a move, so set the inMove flag, and record the
            % index as the start of a move.
            inMove = true;
        end 
    else
        % threshold not exceeded

        if inMove
            inMove = false;
        end
    end

end

figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, smoothedAcc);
title("Smoothed Acc vs time (red), over Position vs time (blue)");

%%

figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, speed);
title("Speed vs time (red), over Position vs time (blue)");

%%
figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, acc);
title("Acc vs time (red), over Position vs time (blue)");

%%
figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, smoothedAcc);
title("Smoothed Acc vs time (red), over Position vs time (blue)");

%%
figure;
plot(DataNorm.time, speed);
yyaxis right;
plot(DataNorm.time, smoothedAcc);
title("Smoothed Acc vs time (red), over Speed vs time (blue)");


%%

for i=1:length(locs)
    xline(locs(i));
end

function TEST1H1 = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  TEST1H1 = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  TEST1H1 = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  TEST1H1 = importfile("P:\NIN Projects\NIN2272\Laser Tracker Data\NIN2272\TEST_1_081223\TEST_1_H.txt", [7, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 13-Dec-2023 13:40:51

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [7, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["x", "y", "z", "time"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
TEST1H1 = readtable(filename, opts);

end