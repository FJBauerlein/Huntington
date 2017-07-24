%% Script to analyze local ER mobility around IBs
% Part 1 of 3

% Felix JB Baeuerlein

% Baeuerlein et al. Cell 2017


%% Input Filenames: 
% [Filename]_1.tif : Reference GFP-image before 10Hz imaging
% [Filename]_2.tif : 10Hz ER-timeseries
% [Filename]_3.tif : Reference GFP-image after 10Hz imaging
% [Filename]_4.tif : 5Hz ER-timeseries (not necessary - only experimental)

%% Output Filenames: 
% Output files are cutouts of the IB with local neighbourhood +-5um around IB

% IB[IB#]_[Filename]_2.tif_Average_IB_[# of IB in FOV].tif : 
%    Time series average of red ER-channel combined with green IB-GFP channel
% IB[IB#]_[Filename]_2.tif_GFP_[# of IB in FOV].tif : 
%    green IB-GFP channel
% IB[IB#]_[Filename]_2.tif_RadialAverage_IB_[# of IB in FOV].tif : 
%    Radial average of the ER and the variance of the ER signal around the IB
% IB[IB#]_[Filename]_2.tif_Timeseries_BC_IB_[# of IB in FOV].tif : 
%    Full time series of ER-channel - after bleach correction
% IB[IB#]_[Filename]_2.tif_Timeseries_IB_[# of IB in FOV].tif : 
%    Full time series of ER-channel - before bleach correction
% IB[IB#]_[Filename]_2.tif_Variance_IB_[# of IB in FOV].tif : 
%    Variance of the full time series of the ER-channel

%% script body

clear all
close all
tic

[FileName,PathName] = uigetfile('*.tif','Select one file in the correct folder');

cd(PathName)
D=dir;
mkdir('Analysis')
IB_count = 0;
fid = fopen('IB_FIJI_Macro.txt','wt');


%% load TimeLapse of 10Hz Imaging
Extension = '2.tif'; % for 10Hz: 2.tif For 5Hz 4.tif 
for n=3:size(D,1)
    tic
    if isempty(strfind(D(n).name, '_4.tif'))==0 && (D(n).isdir==0) % loop checking all time-series
        fname_series = [D(n).name(1:size(D(n).name,2)-5) Extension];
        fname_GFP1 = [D(n).name(1:size(D(n).name,2)-5) '1.tif']; % Reference GFP-image before 10Hz imaging
        info = imfinfo(fname_series);
        Series_10 = [];
        numberOfImages = length(info);
        for k = 1:numberOfImages
            currentImage = imread(fname_series, k, 'Info', info);
            Series_10(:,:,k) = currentImage;
        end
        GFP1.original = imread(fname_GFP1);
        disp(['-- File: ' fname_series ' loaded --'])
        toc


        %% Identification of IB's in the image
        GFP1.IB = GFP1.original.*sign(sign(GFP1.original - 12*median(median(GFP1.original)))); % thresholds the GFP image to detect IB's
        GFP1.props = regionprops(logical(GFP1.IB),'centroid','Area','BoundingBox','Eccentricity','MajorAxisLength', 'MinorAxisLength');  % measures location and properties of IB's       
        disp(['- ' num2str(size(GFP1.props,1)) ' IBs found in image'])
        center = round(cat(1,GFP1.props.Centroid)); % center coordinates of all IB's in the image

        %% calculate per IB 
        xSize = 70; % with OPS 102nm = +-5um around IB
        ySize = 70; % with OPS 102nm = +-5um around IB
        Length = size(Series_10,3); % Length of time series
        
        for m = 1:size(center,1)
            IB.original(:,:,:,m) = Series_10(center(m,2)-xSize:center(m,2)+xSize,center(m,1)-ySize:center(m,1)+ySize,:); % cut out IB's from ER series
            IB_count = IB_count + 1;


            %% calculate average of ER time series and GFP
            GFP_shift_y = 0;
            GFP_shift_x = 0;
            IB.average(:,:,1,m) = mean(IB.original(:,:,:,m),3); % red channel - ERmCherry
            IB.average(:,:,2,m) = GFP1.original(center(m,2)-xSize+GFP_shift_x:center(m,2)+xSize+GFP_shift_x,center(m,1)-ySize+GFP_shift_y:center(m,1)+ySize+GFP_shift_y,:); % green channel - GFP
            % IB.average(:,:,3,m) blue channel in next section
            GFP1.IBcut(:,:,m) = IB.average(:,:,2,m); % cut out IB's from GFP image


            %% define IB boarder
            GFP_decay = 0.4; % Threshold of max GFP-Fluorescence where the IB boarder is defined
            GFP1.IB_Boarder(:,:,m) = sign((GFP1.IBcut(:,:,m) - GFP_decay*max(max(GFP1.IBcut(:,:,m))))); % defines the boarder pixels of the IB
            GFP1.IB_Boarder(GFP1.IB_Boarder<0) = 0;
            GFP1.IB_Boarder(:,:,m) = bwperim(GFP1.IB_Boarder(:,:,m),4);
            IB.average(:,:,3,m) = GFP1.IB_Boarder(:,:,m); % blue channel - IB boarder


            %% ---------   Calculate the ER-Mobility -----------------
            % will be calculated in part2
            
            
            %% save files
            for k = 1:size(IB.original,3)
                imwrite(uint16(IB.original(:,:,k,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_Timeseries_IB_' num2str(m) '.tif'],'writemode','append')
            end
            imwrite(uint16(IB.average(:,:,:,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_Average_IB_' num2str(m) '.tif'])
            imwrite(uint16(GFP1.IBcut(:,:,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_GFP_IB_' num2str(m) '.tif'])

            toc
            disp(['- files written'])
            
            %% write FIJI-Macro to BleachCorrect there
fprintf(fid, ['open("' strrep(PathName,'\','\\\\') 'Analysis\\\\IB' num2str(IB_count,'%03i') '_' [D(n).name(1:size(D(n).name,2)-5) Extension] '_Timeseries_IB_' num2str(m) '.tif");\n']);
fprintf(fid, ['selectWindow("IB' num2str(IB_count,'%03i') '_' [D(n).name(1:size(D(n).name,2)-5) Extension] '_Timeseries_IB_' num2str(m) '.tif");\n']);
fprintf(fid, ['run("Bleach Correction", "correction=[Exponential Fit]");']);
fprintf(fid, ['selectWindow("y = a*exp(-bx) + c");\n']);
fprintf(fid, ['close();\n']);
fprintf(fid, ['selectWindow("IB' num2str(IB_count,'%03i') '_' [D(n).name(1:size(D(n).name,2)-5) Extension] '_Timeseries_IB_' num2str(m) '.tif");\n']);
fprintf(fid, ['close();\n']);
fprintf(fid, ['selectWindow("DUP_IB' num2str(IB_count,'%03i') '_' [D(n).name(1:size(D(n).name,2)-5) Extension] '_Timeseries_IB_' num2str(m) '.tif");\n']);
fprintf(fid, ['saveAs("Tiff", "' strrep(PathName,'\','\\\\') 'Analysis\\\\IB' num2str(IB_count,'%03i') '_' [D(n).name(1:size(D(n).name,2)-5) Extension] '_Timeseries_BC_IB_' num2str(m) '.tif");\n']);
fprintf(fid, ['close();\n']);
         
        end
    end
end

fclose(fid)
toc

% run FIJI-Macro before starting part 2 of the ER_mbilitly_analysis script!