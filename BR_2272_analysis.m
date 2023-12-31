% 
%  Laser tracker import and processing script
%  Author: Ben Rae
%  Last edited: 13/12/2023
%
%  Notes: Currently, the script can process the data, calculate
%  derivatives, and plot in stackedplots (good way to correlate data with
%  movements)
%
%  I wanted to try and automatically split moves based on measuring when
%  acceleration or jerk exceed a threshold - however doing this
%  automatically is a bit tricky. 
%
%  Instead I think it might make sense to use thresholding to detect
%  changing jerk, i.e. beginning of any meaningful data for any given
%  machine move, then separate data manually into individual moves (if we
%  want to split them at all). 
% 
%  Using the jerk for each axis individually might not be much use - next
%  I'm wondering whether to calculate displacement/speed/acceleration/jerk
%  in 3D
%

clc
clear
close all

% Import data using automated importer
Data = importfile('P:\NIN Projects\NIN2272\Laser Tracker Data\NIN2272\TEST_1_081223\TEST_1_H.txt');
%Data = importfile('A:\namrc\2272_analysis\Test data\TEST_3_H.txt');

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

%% Calculate derivatives
                  
X_Spd = smoothdata(gradient(Data.x),"movmean",200);   % 1st derivative / speed for X             
Y_Spd = smoothdata(gradient(Data.y),"movmean",200);   % 1st derivative / speed for Y
Z_Spd = smoothdata(gradient(Data.z),"movmean",200);   % 1st derivative / speed for Z

X_Acc = smoothdata(gradient(X_Spd),"movmean",400);    % 2nd derivative / acceleration for X
Y_Acc = smoothdata(gradient(Y_Spd),"movmean",400);    % 2nd derivative / acceleration for Y
Z_Acc = smoothdata(gradient(Z_Spd),"movmean",400);    % 2nd derivative / acceleration for Z

X_Jrk = smoothdata(gradient(X_Acc),"gaussian", 1000); % 3rd derivative / jerk for X
Y_Jrk = smoothdata(gradient(Y_Acc),"gaussian", 1000); % 3rd derivative / jerk for Y
Z_Jrk = smoothdata(gradient(Z_Acc),"gaussian", 1000); % 3rd derivative / jerk for Z

%% Convert data table into timetable (allows extra processing)

outfmt = 'hh:mm:ss.SSS'; 
D = duration(0, 0, Data.time,'Format',outfmt); % convert time into duration format using above format

% Create overall data table with all data
T = table(D,Data.x,Data.y,Data.z,X_Spd,X_Acc,X_Jrk,Y_Spd,Y_Acc,Y_Jrk,Z_Spd,Z_Acc,Z_Jrk,'VariableNames',["Time","X_Pos","Y_Pos","Z_Pos","X_Spd","X_Acc","X_Jrk","Y_Spd","Y_Acc","Y_Jrk","Z_Spd","Z_Acc","Z_Jrk"]);
DataTable = table2timetable(T); 

% Create separated data tables for individual axes
X_Data = table(D, Data.x, X_Spd, X_Acc, X_Jrk, 'VariableNames',["Time", "Position", "Speed", "Acceleration", "Jerk"]);
X_Data = table2timetable(X_Data);
Y_Data = table(D, Data.y, Y_Spd, Y_Acc, Y_Jrk, 'VariableNames',["Time", "Position", "Speed", "Acceleration", "Jerk"]);
Y_Data = table2timetable(Y_Data);
Z_Data = table(D, Data.z, Z_Spd, Z_Acc, Z_Jrk, 'VariableNames',["Time", "Position", "Speed", "Acceleration", "Jerk"]);
Z_Data = table2timetable(Z_Data);

%% Graphing of overall movements (3D + separated X, Y, Z)

t = tiledlayout(3,6);
ax1 = nexttile(1, [3 3]);
hold on;
f3d = plot3(DataTable.X_Pos,DataTable.Y_Pos,DataTable.Z_Pos);
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

ax2 = nexttile(4, [3 3]);
s = stackedplot(DataTable, ["X_Pos", "Y_Pos", "Z_Pos"]);
title("Individual axis positions");
grid on;

%% Graphing of X/Y/Z movement components
figure;
t1 = tiledlayout(3,6);

ax3 = nexttile(1, [3 2]);
s = stackedplot(DataTable, ["X_Pos", "X_Spd", "X_Acc", "X_Jrk"]);
title("Data for linear move in X axis, high feedrate (4000 mm/min)");
grid on;

ax4 = nexttile(3, [3 2]);
s = stackedplot(DataTable, ["Y_Pos", "Y_Spd", "Y_Acc", "Y_Jrk"]);
title("Data for linear move in Y axis, high feedrate (4000 mm/min)");
grid on;

ax4 = nexttile(5, [3 2]);
s = stackedplot(DataTable, ["Z_Pos", "Z_Spd", "Z_Acc", "Z_Jrk"]);
title("Data for linear move in Y axis, high feedrate (4000 mm/min)");
grid on;

%% Graphing of X/Y/Z movement components but overlaid on same axes
figure;
s = stackedplot(X_Data,Y_Data,Z_Data,"Title","Linear move components, high feedrate (4000 mm/min)");
grid on;

%% Calculate 3D displacement (speed)
figure;
disp = zeros(length(DataTable.Time),1);
for i=2:length(DataTable.Time)
    P1 = [X_Data.Position(i-1) Y_Data.Position(i-1) Z_Data.Position(i-1)];
    P2 = [X_Data.Position(i) Y_Data.Position(i) Z_Data.Position(i)];
    disp(i) = norm(P2 - P1);
end

3d

%% Try to separate the movements based on jerk. When jerk is not zero, mark cut in data.
% Currently this won't work

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
figure;
s = stackedplot(DataTable, ["X_Pos", "X_Spd", "X_Acc", "X_Jrk"]);
title("Data for linear move in X axis, high feedrate (4000 mm/min)");

%% Stacked plot to show multiple data
figure;
s = stackedplot(DataTable, ["X_Acc", "Y_Acc", "Z_Acc"]);
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