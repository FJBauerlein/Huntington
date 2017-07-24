%% script to transform data and calculate the Persistance Length of filaments from Amira 6.0.1

% Data needed: Export from Amira 6.0.1 - attribute file in .xml
% open in Excel and save as .xlsx

% Felix JB Baeuerlein

% Baeuerlein et al. Cell 2017


%% script body

% Load data
clear all
clc
[filename, pathname] = uigetfile({'*.xlsx'},'File Selector');
disp(['----- Processing Tomogram: ' filename ' -----'])

% Transform dataformat from Amira 6.0.1 in 5.4.5 format
tic
disp('1) Starting Amira Data transformation ...')
Amira_lineset_conversion
toc
disp('Data transformed and saved!')

% copy files into folder
f = regexp(filename, '.xlsx', 'split');
fname = f{1,1}; 
mkdir(fname);
movefile(filename,fname);
movefile(([filename '_Points.csv']),fname);
cd(fname);
disp('Data moved to Folder!')


% calculate motion direction
tic
disp('2) Starting PCA for general motion direction ...')
filaments=Points; 
fil_motion_direction_m2;
% movefile([filename '_with_motion_direction.csv'],fname);
toc
disp('Motion Direction calculated!')


% interpolate equidistant points on each filament
tic
disp('3) Starting equidistant interpolation ...')
searchdist = 50;
fil_analysis_nofigs_m2;
toc
disp('Equidistant interpolation done!')


% calculate the persistance length
tic
disp('4) Starting Persistence Length calculation ...')
Persistencelength_160712
toc
disp('Persistance Length calculated!')

% save final coordinate file
save([fname '_final_Points.mat'],'filaments');

disp('Final points list saved!')


% check Algorithm funcionality
check=filaments(2:end,1)-filaments(1:end-1,1)-1;
disp(['Number of lost points: -- ' num2str(sum(check)) ' -- from total ' num2str(max(filaments(:,1))) ' points.'])

% visualizes the first 50 filaments
% figure(99);
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