%% Script to analyze local ER mobility around IBs
% Part 2 of 3

% Felix JB Baeuerlein

% Baeuerlein et al. Cell 2017


%% script body

close all

[FileName,PathName] = uigetfile('*.tif','Select one file in the correct folder');

cd(PathName)
D=dir;
IB_count = 0;

tic
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
            IB_count = IB_count + 1;
            fname_BleachCorr = (['Analysis/IB' num2str(IB_count,'%03i') '_' D(n).name(1:size(D(n).name,2)-5) Extension '_Timeseries_BC_IB_' num2str(m) '.tif']);
            info = imfinfo(fname_BleachCorr);
            IB.original = [];
            numberOfImages = length(info);
            for k = 1:numberOfImages
                currentImage = imread(fname_BleachCorr, k, 'Info', info);
                IB.original(:,:,k) = currentImage;
            end
            

            %% calculate average of ER time series and GFP
            GFP_shift_y = 0;
            GFP_shift_x = 0;
            IB.average(:,:,1,m) = mean(IB.original(:,:,:),3); % red channel - ERmCherry
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
            tic
            %% Strategy I: Simple Variance of time-series
            IB.Var(:,:,m) = nanvar(IB.original(:,:,:),1,3); % calculates the variance of the time-series
            IB.Var_norm(:,:,m) = IB.Var(:,:,m)./IB.average(:,:,1,m); % pixelwise normalization to mean Intensity

            %% Strategy II: Variance of the Difference between two images spaced by 'stepsize' (w or w/o normalization)
%             stepsize = 4; % distance between images in stack
%             IB.Delta(:,:,:,m) = IB.original(:,:,1:stepsize:(Length-stepsize)) - IB.original(:,:,2:stepsize:Length-(stepsize-1)); % calculates the difference of the images
%             % IB.Delta(:,:,:,m) = IB.original(:,:,1:Length-stepsize,m) - IB.original(:,:,stepsize+1:Length,m); 
%             % IB.Delta(:,:,:,m) = IB.Delta(:,:,:,m) ./ max(IB.original(:,:,1:stepsize:(Length-stepsize),m),IB.original(:,:,2:stepsize:Length-(stepsize-1),m)); % normalizes the Delta to the pixels max intensity
%             IB.Var(:,:,m) = nanvar(IB.Delta(:,:,:,m),1,3); % calculates the variance of the the image-difference


            %% Strategy III: CrossCorrelation of subimages over time
%             subimg_size = 9;
%             px = floor(subimg_size/2);
%             for i=2:size(IB.original,3) % loop for timepoints
%                 CC(i) = max(max(corr2(IB.original(:,:,i-1),IB.original(:,:,i)))); % a) 1st vs ith  b) i-1 vs i
%             end
%             plot(1:size(CC,2),CC)
%             
%             tic
%             IB.CC = NaN(size(IB.original,1),size(IB.original,2),size(IB.original,3));
% %             IB.original=gpuArray(IB.original);
%             for p=px+1:size(IB.original,1)-px % loop for x
%                 for q=px+1:size(IB.original,2)-px % loop for y
%                     for i=2:size(IB.original,3) % loop for timepoints
%                         IB.CC(p,q,i) = max(max(corr2(IB.original(p-px:p+px,q-px:q+px,i-1),IB.original(p-px:p+px,q-px:q+px,i)))); % a) 1st vs ith  b) i-1 vs i
%                     end
%                 end
%             end
%             IB.Var(:,:,m) = nanvar(IB.CC,1,3);
%             imwrite(uint16(IB.Var(:,:,m)*1e5),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_VarianceCC_IB_' num2str(m) '.tif'])
% %             imshow(IB.Var(:,:,m), [0 0.025])
%             toc
%             


            %% save and visualize Variance-Result
            
            imwrite(uint16(IB.average(:,:,:,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_Average_IB_' num2str(m) '.tif'])
            imwrite(uint16(IB.Var(:,:,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_Variance_IB_' num2str(m) '.tif'])
            imwrite(uint16(IB.Var_norm(:,:,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_VarianceNormPx_IB_' num2str(m) '.tif'])
            imwrite(uint16(GFP1.IBcut(:,:,m)),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_GFP_IB_' num2str(m) '.tif'])

            toc
            disp(['- ER-Mobility calcuated'])


            %% Calculate rotational average
            % presetting the polar coordinate transformation
            [xs ys] = size(IB.Var(:,:,m));
            f2 = -xs/2:xs/2-1;
            f1 = -ys/2:ys/2-1;
            [XX YY] = meshgrid(f1,f2);
            [t r] = cart2pol(XX,YY);
            t = round(t./pi.*180); t(t<0)=t(t<0)+360; % change into degrees
            if mod(xs,2)==1 || mod(ys,2)==1
                r = round(r)-1;
            else
                r = round(r);
            end

            %% define non-cytosolic regions (Ncl, extracellular)
            excl_threshold = 60; % Strategy I : 50-60 - Strategy II : 90-100 - Strategy III : 0
            Var = IB.Var(:,:,m);
            ER = IB.average(:,:,1,m);
            GFP =  GFP1.IBcut(:,:,m);
            ER(Var<(excl_threshold-(GFP>0.1*max(max(GFP1.IBcut(:,:,m))))*excl_threshold)) = NaN; % sets all pixels except the IB to NaN that are in the Ncl or extracellular
            Var(Var<(excl_threshold-(GFP>0.1*max(max(GFP1.IBcut(:,:,m))))*excl_threshold)) = NaN; % sets all pixels except the IB to NaN that are in the Ncl or extracellular
            Var_normPx = Var./ER;
%             imwrite(logical(Var(Var>0)==1),['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_VarianceMask_IB_' num2str(m) '.tif'])
            Boarder = IB.average(:,:,3,m);
            
            for angle = 1:1:360
                if angle == 1
                    last_dist = r.*Boarder;
                    last_dist(last_dist==0)=NaN;
                    last_dist = round(nanmean(nanmean(last_dist)));
                end
                if max(max(r(t==angle).*Boarder(t==angle))) > 0
                    boarder_dist(angle) = max(max(r(t==angle).*Boarder(t==angle)));
                    last_dist = max(max(r(t==angle).*Boarder(t==angle)));
                else
                    boarder_dist(angle) = last_dist;
                end
                r(t==angle) = r(t==angle)-boarder_dist(angle);
            end

            %% average rotationally
            R_neg = 20; % for claculating values inside the IB - 1unit = 100nm
            for radius = 1:floor((min(xs,ys)+R_neg)/2) % min(min(r)):max(max(r))
                Radial.Average_Var(radius,IB_count) = nanmean(Var(r==radius-R_neg));
                Radial.Average_ER(radius,IB_count) = nanmean(ER(r==radius-R_neg));
                Radial.Average_Var_norm(radius,IB_count) = Radial.Average_Var(radius,IB_count)./Radial.Average_ER(radius,IB_count);
                Radial.Average_Var_normPx(radius,IB_count) = nanmean(Var_normPx(r==radius-R_neg));
                Radial.Average_GFP(radius,IB_count) = nanmean(GFP(r==radius-R_neg));
                Radial.Radius(radius) = (radius-R_neg)*0.1;
            end

            figure;
            plot(Radial.Radius(:),Radial.Average_Var(:,IB_count),'b'); 
            hold on
            plot(Radial.Radius(:),Radial.Average_ER(:,IB_count),'r');
            plot(Radial.Radius(:),Radial.Average_GFP(:,IB_count),'g');
            plot(Radial.Radius(:),Radial.Average_Var_norm(:,IB_count)*nanmean(nanmean(Radial.Average_ER)),'c')
            plot(Radial.Radius(:),Radial.Average_Var_normPx(:,IB_count)*nanmean(nanmean(Radial.Average_ER)),'c--')
            ylim([0 1.5*max(max(Radial.Average_ER(:,IB_count)))])
            xlabel('Distance in um');
            ylabel('Variance');
            saveas(gcf,['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_RadialAverage_IB_' num2str(m) '.png'])
            hold off
        end
    end
end

save('RadialAverages.mat','Radial')
toc