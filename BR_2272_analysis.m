% 
%  Laser tracker import and processing script
%  Author: Ben Rae
%  Last edited: 13/12/2023
%

clc
clear
close all

% Import data using automated importer
Data = importfile('P:\NIN Projects\NIN2272\Laser Tracker Data\NIN2272\TEST_1_081223\TEST_1_H.txt');

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
% DataNorm.x = abs(DataNorm.x);
% DataNorm.y = abs(DataNorm.y);
% DataNorm.z = abs(DataNorm.z);

DataNorm.x = DataNorm.x;
DataNorm.y = DataNorm.y;
DataNorm.z = DataNorm.z;

%% Graphing of overall movements (3D + separated X, Y, Z)

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

%% Calculate derivatives

speed = gradient(DataNorm.x);                       % 1st derivative / speed
smoothedSpeed = smoothdata(speed,"movmean",200);    % smoothing speed data
smoothedSpeedAbs = abs(smoothedSpeed);              % absolute set of speed data (all data positive)
acc = gradient(speed);                              % 2nd derivative / acceleration 
smoothedAcc = smoothdata(acc,"movmean",400);        % smoothing acc data
smoothedAccAbs = abs(smoothedAcc);                  % absolute set of acc data (all data positive)
jerk = gradient(acc);                               % 3rd derivative / jerk
smoothedJerk = smoothdata(jerk, "gaussian", 1000);  % smoothing jerk data

%% Try to separate the movements based on speed. When speed is not zero, mark cut in data.

threshold = 0.0007; % difference from zero that implies acceleration/movement
inMove = false; % flag indicating whether we're in a move, i.e. not static
startIndices = 0;
startIndicesCounter = 1;
endIndices = 0;
endIndicesCounter = 1;

for i=1:length(smoothedSpeedAbs) 
    if smoothedSpeedAbs(i) > threshold
        % postive threshold exceeded
        if inMove
            % we're already in a move, do nothing.            
        else 
            % we're not in a move, so set the inMove flag, and record the
            % index as the start of a move.
            inMove = true;
            startIndices(startIndicesCounter) = i;
            startIndicesCounter = startIndicesCounter + 1;
        end        
    else
        % threshold not exceeded
        if inMove
            % we're in a move, and threshold is now not exceeded - move is
            % concluded, so record end index.
            inMove = false;
            endIndices(endIndicesCounter) = i;
            endIndicesCounter = endIndicesCounter + 1;
        else
            % we're not in a move but threshold is not exceeded, do
            % nothing.
        end
    end
end

disp("Recorded " + length(startIndices) + " move start points, and " + length(endIndices) + " move end points.");

figure;
plot(DataNorm.time, DataNorm.x);
% yyaxis right;
% plot(DataNorm.time, smoothedSpeedAbs);
title("Position data with start/end points of moves indicated");

for i=1:length(startIndices)
    xline(DataNorm.time(startIndices(i)),'g');
end

for i=1:length(endIndices)
    xline(DataNorm.time(endIndices(i)),'r');
end

%% Stacked plot to show multiple data
outfmt = 'hh:mm:ss.SSS';
D = duration(0, 0, DataNorm.time,'Format',outfmt);
T = table(D,DataNorm.x,smoothedSpeed,smoothedAcc,smoothedJerk,'VariableNames',["Time","X Position","Speed","Acceleration","Jerk"]);
XData = table2timetable(T);
figure;
s = stackedplot(XData);
title("Data for linear move in X axis, high feedrate (4000 mm/min)");

%% Plot overlay lines showing moves
ax = findobj(s.NodeChildren, 'Type', 'Axes');
for j=1:length(ax)
    for i=1:length(startIndices)
        xline(ax(j),XData.Time(startIndices(i)),'g');
    end
    
    for i=1:length(endIndices)
        xline(ax(j),XData.Time(endIndices(i)),'r');
    end
end

%% Below - individual sections for plotting data, not required to run! 
% Just for testing/visualising

%% Plot speed over position
figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, speed);
title("Speed vs time (red), over Position vs time (blue)");

%% Plot smoothedSpeedAbs over position
figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, smoothedSpeedAbs);
title("Smoothed, Absolute Speed vs time (red), over Position vs time (blue)");

%% Plot acc over position
figure;
plot(DataNorm.time, DataNorm.x);
yyaxis right;
plot(DataNorm.time, smoothedAcc);
title("Smoothed Acc vs time (red), over Position vs time (blue)");

%% Plot acc over speed
figure;
plot(DataNorm.time, speed);
yyaxis right;
plot(DataNorm.time, smoothedAcc);
title("Smoothed Acc vs time (red), over Speed vs time (blue)");

%% Plot jerk over position
figure;
plot(DataNorm.time, smoothedJerk);
yyaxis right
plot(DataNorm.time, DataNorm.x);
title("Smoothed Jerk over Position");

%% Plot jerk over acc
figure;
plot(DataNorm.time, smoothedJerk);
yyaxis right
plot(DataNorm.time, smoothedAcc);
title("Smoothed Jerk over Acc");

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