function []=RSDA(tile_x, tile_y)

%% function to detect aggregates stained with mCherry-Ub in LM-data from CorrSight

% Felix JB Baeuerlein 
% 17.March 2015

% Baeuerlein et al. Cell 2017



%% Parameters
grd_threshold = 0.1;
grd_vector = 1.0;
ring_area = 50;
agg_area_min = 30;
agg_area_max = 1000;


% --- xlwrite for Mac ---
% %% Initialisation of POI Libs
% % Add Java POI Libs to matlab javapath
% javaaddpath('poi_library/poi-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
% javaaddpath('poi_library/xmlbeans-2.3.0.jar');
% javaaddpath('poi_library/dom4j-1.6.1.jar');
% javaaddpath('poi_library/stax-api-1.0.1.jar');

folder_name=uigetdir
cd(folder_name)

if exist('Originals_SD') == 7
    cd Originals_SD
    if size(dir,1) > 3
        disp('----- Abort: Files already modified! -----')
        cd ..;
        return
    end
    cd ..;
end

%% reading images into memory
disp('-----------------------------------------------------------------')
disp('---------- Aggregate Detecetion Algorithm initiated ...----------')
disp('-----------------------------------------------------------------')
disp('---------- Reading images ...----------')
i=0;
images=struct;
tic
for x=1:tile_x
    for y=1:tile_y
        i=i+1;
%        images(x,y).bf=mat2gray(imread(['Tile_00' num2str(x) '-00' num2str(y) '-000_0.tif']));
        images(x,y).gfp=mat2gray(imread(['Tile_00' num2str(x) '-00' num2str(y) '-000_0.tif']));
        images(x,y).name=['Tile_00' num2str(x) '-00' num2str(y) '-000_0.tif'];
        disp([' Image x:' num2str(x) ' y:' num2str(y) ' (' num2str(i) '/' num2str(tile_x*tile_y) ') loaded... remaining time: ' num2str((toc)/(i)*(tile_x*tile_y-i)) ' s'])
    end
end
toc



%% detection of large/small aggregates
disp('---------- Aggregate Detection started ... ----------')
i=0; j=1;
sizex=1:size(images(1,1).gfp,1); sizex=sizex';
sizey=1:size(images(1,1).gfp,2); sizey=sizey';
ox=ones(1,length(sizex));
oy=ones(1,length(sizey));
tic
for x=1:tile_x
    for y=1:tile_y
        i=i+1;
%% detection of aggregates in fluorescence-images
        image = mat2gray(images(x,y).gfp);   % data conversion
        [m,n] = size(image);
        [X Y] = meshgrid(1:n,1:m);
        [grdx,grdy] = gradient(image,grd_vector);   % gradient in x and y
        grds = abs(grdx)+abs(grdy);   % gradient in 2D

%% select structures in the right size
        grds = grds-grd_threshold;   % selection of high gradient structures
        grds(grds<0) = 0;   % set flat stuff to zero
        grds = bwlabel(grds,4);   % label structures that belong together
        st = regionprops( grds, 'Centroid', 'Area', 'Eccentricity' );   % get parameters
        toosmall = [st(:).Area]<ring_area;   % get rid of too small structures
        for k = size(toosmall,2):-1:1 
            if toosmall(1,k) == 1
                grds(grds==k)=0;   % set gradient of too small structures to zero
                st(k) = [];   % delete entry from list
            end
        end

%% select ring structures
        grds = bwmorph(grds,'skel',inf);
        grds = -grds+1;
        grds(grds<1) = 0;  % flat areas are 1
        grds = bwlabel(grds,4);  % ring structures are 0

        Aggregates(x,y).props = regionprops( grds, 'Centroid', 'Area', 'Boundingbox', 'MajorAxisLength', 'MinorAxisLength' );
        toosmall2 = [Aggregates(x,y).props.Area]<agg_area_min;  % too small areas
        toobig2 = [Aggregates(x,y).props.Area]>agg_area_max;  % background, apoptitic cells
        too2 = toosmall2 + toobig2;  % list for deletion
        for k = size(too2,2):-1:1
            if too2(1,k) == 1
                grds(grds==k)=0;
                Aggregates(x,y).props(k,:)=[];
            end
        end
        Aggregates(x,y).Area=cat(1,Aggregates(x,y).props.Area);


%% write Aggregate List for each image
        Aggregates(x,y).Agg=[];
        for n=1:length(Aggregates(x,y).Area)
                    Aggregates(x,y).Agg(n,1)=x;
                    Aggregates(x,y).Agg(n,2)=y;
        end
        if isempty(Aggregates(x,y).Agg)==0
            Aggregates(x,y).Agg(:,5)=Aggregates(x,y).Area(:,1);
            Aggregates(x,y).Agg(:,7:8)=round(cat(1,Aggregates(x,y).props.Centroid));
            Aggregates(x,y).Agg(:,9:12)=round(cat(1,Aggregates(x,y).props.BoundingBox));
            Aggregates(x,y).Agg(:,4)=sqrt(Aggregates(x,y).Agg(:,5)/pi())*2*0.364;

%% size exclusion of detected signals
%             toosmall2 = [Aggregates(x,y).Area]<20;  % too small areas
%             toobig2 = [Aggregates(x,y).Area]>1000;  % background, apoptitic cells
%             too2 = toosmall2 + toobig2;  % list for deletion
%             for k = size(too2,2):-1:1
%                 if too2(1,k) == 1
%                     grds(grds==k)=0;
%                     Aggregates(x,y).Agg(k,:)=[];
%                 end
%             end
%             AggList_toosmall=Aggregates(x,y).Agg(:,5)<20; % which Detections are too small for Aggregates
%             Aggregates(x,y).Agg(AggList_toosmall,:)=[]; % delete too small detections
%             AggList_toolarge=Aggregates(x,y).Agg(:,5)>5000; % which Detections are too large for Aggregates
%             Aggregates(x,y).Agg(AggList_toolarge,:)=[]; % delete too large detections
        end
        AggList.value=cat(1,Aggregates.Agg); % writes a cumulated list of all aggregates

%% display progress
        disp([' Aggregate Detection in x:' num2str(x) ' y:' num2str(y) ' (' num2str(i) '/' num2str(tile_x*tile_y) ') remaining time: ' num2str((toc)/(i)*(tile_x*tile_y-i)) ' s'])
%         images(x,y).bf=1; images(x,y).mask=1;
    end
end
AggList.name={'x Tile' 'y Tile' 'Score' 'Diameter [?m]' 'Area [px]' 'Distance to Gridbar [?m]' 'x Position' 'y Position' 'x BBox LU' 'y BBox LU' 'x BBox Length' 'y BBox Length'};



%% size exclusion of detected signals
AggList_final=AggList.value;
% AggList_toosmall=AggList_reduced(:,4)<100; % which Detections are too small for Aggregates
% AggList_reduced(AggList_toosmall,:)=[]; % delete too small detections
% AggList_toolarge=AggList_reduced(:,4)>5000; % which Detections are too large for Aggregates
% AggList_reduced(AggList_toolarge,:)=[]; % delete too large detections
AggList_final=sortrows(AggList_final,-5);
toc


%% write variables to workspace
assignin('base','images', images)
assignin('base','Aggregates', Aggregates)
assignin('base','AggList', AggList)
assignin('base','AggList_final', AggList_final) 

%% make Excel-File with table and variable names
% AggList_final(:,3)=round(AggList_final(:,3)*100)/100;
AggList_final(:,4)=round(AggList_final(:,4)*10)/10;
% AggList_final(:,6)=round(AggList_final(:,6)*10)/10;
Table=[AggList.name;num2cell(AggList_final)];
% writetable(cell2table(Table),['Aggregate_List('  folder_name(max(findstr(folder_name,'/'))+1:end) ').csv']);
xlwrite(['Aggregate_List('  folder_name(max(findstr(folder_name,'/'))+1:end) ')' ],Table); % xlswrite for Mac
% xlswrite(['Aggregate_List('folder_name(max(findstr(folder_name,'/'))+1:end) ')' ],Table); % doesn't work with Mac


%% mark aggregates, move originals and save modified GFP-images
disp('---------- Modification of fluorescence images started ... -------------')
i=0;
tic
mkdir Originals_SD
mkdir modif
boxsize=100/2;
for x=1:tile_x
    for y=1:tile_y
        if isempty(Aggregates(x,y).Agg)==0
            images(x,y).gfp=imread(['Tile_00' num2str(x) '-00' num2str(y) '-000_0.tif']);
            info=imfinfo(['Tile_00' num2str(x) '-00' num2str(y) '-000_0.tif']);
            intensity=max(max(images(x,y).gfp));
            for n=1:length(Aggregates(x,y).Agg(:,1))
                i=i+1;
% mark aggregates in fluorescence-images
                y1=round(Aggregates(x,y).Agg(n,9))-boxsize;
                if y1<0
                    y1=1;
                end
                y2=round(Aggregates(x,y).Agg(n,9))+round(Aggregates(x,y).Agg(n,11))+boxsize;
                if y2>1024;
                    y2=1024;
                end
                x1=round(Aggregates(x,y).Agg(n,10))-boxsize;
                if x1<0
                    x1=1;
                end
                x2=round(Aggregates(x,y).Agg(n,10))+round(Aggregates(x,y).Agg(n,12))+boxsize;
                if x2>1344;
                    x2=1344;
                end
                images(x,y).gfp(x1:x2,y1:y1+1)=intensity*100;
                images(x,y).gfp(x1:x2,y2-1:y2)=intensity*100;
                images(x,y).gfp(x1:x1+1,y1:y2)=intensity*100;
                images(x,y).gfp(x2-1:x2,y1:y2)=intensity*100;
% save images to Maps-Folder, copy originals in folder originals_SD
            end
            movefile(images(x,y).name,'Originals_SD/') % moves all Tile-Files in the 'original' Folder
            imwrite(images(x,y).gfp, ['modif/' images(x,y).name],'tif','Compression','none','RowsPerStrip',1344,'Resolution',96,'Description',info.ImageDescription)
            imwrite(images(x,y).gfp, images(x,y).name,'tif','Compression','none','RowsPerStrip',1344,'Resolution',96,'Description',info.ImageDescription)
            disp([' Aggregates marked in Tile x:' num2str(x) ' y:' num2str(y) ' (' num2str(i) '/' num2str(length(AggList.value)) ') remaining time: ' num2str((toc)/(i)*(length(AggList.value)-i)) ' s'])
        end
    end
end
toc

disp(['---------- N = ' num2str(size(AggList_final,1)) ' Aggregates found ! ---------------------------------'])
disp('---------- Original fluorescence images moved to folder "Originals_SD" -----')
disp('---------- Excel table of Aggregates saved. ----------------------------')
