%% Amira_lineset_conversion


% Felix JB Baeuerlein

% Baeuerlein et al. Cell 2017

%% script body

% procedure to convert new Amira 6.0.1 Filament Tracing Data into the 
% Amira 5.4.0 format:
% filaments (:,1) = point index (from 0 to the last index, 1 index per point, independently of the line number)
% filaments (:,2) = coord X of the point
% filaments (:,3) = coord Y of the point
% filaments (:,4) = coord Z of the point
% filaments (:,5) = line index (start with 0)
% filaments (:,6) = point index along a line (restart to 0 for each new line)
% filaments (:,7) = CC
% filaments (:,8) = Sim
% filaments (:,9) = theta (from the Y-axis to the X-axis, range: 0 .. 180)
% filaments (:,10) = psi (angle from the XY plane to the Z, range: 0..90)
% filaments (:,11) = line length

% Load data
[filename, pathname] = uigetfile({'*.xlsx'},'File Selector');
% filename='t79_C70_L40_S50.attributegraph.xlsx';
[~, ~, Seg] = xlsread([pathname filename],'Segments');
[~, ~, Points] = xlsread([pathname filename],'Points');

N = 2; % counter for Segment-file line
L = 10; % if Amiradata in ? - to convert to nm


for i= 2:length(Points)
    
    if N < length(Seg)-1
            if min(str2num(Seg{N+1,10})) == i-2
                N=N+1;
            end
            else
                N=length(Seg)-1;
    end
    
    if N < length(Seg)-1
            if Points{i,1} < min(str2num(Seg{N+1,10}))
                Points{i,2} = Points{i,2}./L; % conversion of Lengthes defined in Amira to nm
                Points{i,3} = Points{i,3}./L; % conversion of Lengthes defined in Amira to nm
                Points{i,4} = Points{i,4}./L; % conversion of Lengthes defined in Amira to nm
                Points{i,5}=Seg{N,7}; % line index (start with 0)
                Points{i,6}=Points{i,1}-(min(str2num(Seg{N,10}))-1); % point index
                Points{i,7}=0;
                Points{i,8}=0;
                Points{i,9}=Seg{N,5}; % theta (from the Y-axis to the X-axis, range: 0 .. 180)
                Points{i,10}=Seg{N,6}; % psi (angle from the XY plane to the Z, range: 0..90)
                Points{i,11}=Seg{N,2}./L; % line length
            end
    else
                Points{i,2} = Points{i,2}/L; % conversion of Lengthes defined in Amira to nm
                Points{i,3} = Points{i,3}/L; % conversion of Lengthes defined in Amira to nm
                Points{i,4} = Points{i,4}/L; % conversion of Lengthes defined in Amira to nm
                Points{i,5}=Seg{N,7}; % line index (start with 0)
                Points{i,6}=Points{i,1}-(min(str2num(Seg{N,10}))-1); % point index
                Points{i,7}=0;
                Points{i,8}=0;
                Points{i,9}=Seg{N,5}; % theta (from the Y-axis to the X-axis, range: 0 .. 180)
                Points{i,10}=Seg{N,6}; % psi (angle from the XY plane to the Z, range: 0..90)
                Points{i,11}=Seg{N,2}./L; % line length
    end
    
end


Points=cell2mat(Points(2:length(Points),:));
csvwrite([filename '_Points.csv'], Points)

Line=cell2mat(Seg(2:end,1:9));
Line2(:,1)=Line(:,1);
Line2(:,2)=Line(:,2)/10;
Line2(:,3:4)=Line(:,5:6);
Line2(:,5)=ones(length(Line2),1)*85;
csvwrite([filename '_Lines.csv'], Line2)


% % visualizes the first 50 filaments
% for ff=0:50
%     ind = find(Points(:,5) == ff); 
%     plot3(Points(ind,2), Points(ind,3), Points(ind,4), '*-');
%     axis equal
%     view([-35.5 14]);
%     box('on');
%     grid('on');
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%   %  legend(ind, '# Filament')
%     pause(0.2)
%     hold on
%     ind = find(filaments(:,5) == ff); 
%     plot3(filaments(ind,2), filaments(ind,3), filaments(ind,4), 'or');
%     pause(0.5)
%     hold off
% end