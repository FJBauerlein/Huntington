%% Script to determine the Persistence length of filaments

% Felix JB Baeuerlein

% Baeuerlein et al. Cell 2017


%% script body


% [filename, pathname] = uigetfile({'*.csv'},'File Selector');
% filaments = csvread([pathname filename]);

stepsize = 1;
samplingdistance = 5.0; % distance between two sampling point in nm
s = samplingdistance*stepsize; % Spacing between sampling points in nm

fil_no = filaments(end,5)+1;

Filament=0;
j = 1;             
for f=0:max(filaments(:,5))                         % for all filaments (:,5) = line index
    ind = filaments(find(filaments(:,5) == f),1);   % find all points corresponding to the same filament
    ind = ind +1;                                   % starts with 0 from amira
    temp = filaments(ind(1):ind(end),:);
        if temp(1,11)>25                           % discard filaments shorter than 34 nm (10 data points)   
              Filament = Filament+1; 

              maxCC = round((max(temp(:,1))-min(temp(:,1)))/2)+min(temp(:,1));               % find middle point of filament
              i_maxCC = find(temp(:,1) == maxCC);             % find index point with maximum CC 
              i_maxCC = i_maxCC(1);                 % go in both directions along the filament

              v1 = temp(i_maxCC,14:16);
              data(Filament,1)= f;
              data(Filament,2)= temp(1,11);
              first = i_maxCC;                            
              last1 = size(temp,1);
              last2 = 1;
              ii = 1;
              for i=first:last1                      %  for each filament: can be adjusted to exclude points at edges. all points: i=1:size(temp,1) 
                    PointID = ii;
                    CorAngles(j,PointID) = subspace (v1',temp(i,14:16)') * (2/pi*90); % Angle between tangents of first and n-th samplingpoints with distance s
                    CorCos(j,PointID) = cos (CorAngles(j,PointID)/(2/pi*90));
                    i=i+1;
                    ii=ii+1;
              end 
              ii = first-last2+1;
              j = j+1;
              for i=last2:first                      %  for each filament: can be adjusted to exclude points at edges. all points: i=1:size(temp,1) 
                    PointID = ii;
                    CorAngles(j,PointID) = subspace (v1',temp(i,14:16)') * (2/pi*90); % Angle between tangents of first and n-th samplingpoints with distance s
                    CorCos(j,PointID) = cos (CorAngles(j,PointID)/(2/pi*90));
                    i=i+1;
                    ii=ii-1;
              end
              j=j+1;
              data(Filament,3)= (first-last2)*s;
              data(Filament,4)= (last1-first)*s;
        else    
        f=f+1;    
        end
end

fil_no = size(data,1);
mean_length = mean(data(:,2));
std_length = std(data(:,2));
max_length = max(data(:,2));
edges=0:10:600;
figure; 
h = histogram(data(:,2),edges);
xlabel('Filament length [nm]')
grid on
hold on

% find 95% Percentile to set upper Length cutoff for fitting
h_cum(1,1) = h.Values(1,1);
for i = 2:length(h.Values)
    h_cum(i,1)=h_cum(i-1,1)+h.Values(1,i); % calculates the cumulative histogram
end
h_cum=h_cum/max(h_cum)*100;
plot(edges(1:end-1), h_cum, 'r-')
P90 = edges(1,find(h_cum(:,1) == max(h_cum(h_cum<90))));
plot(edges(1,find(h_cum(:,1) == max(h_cum(h_cum<90)))) , 90, 'ro')
hold off

saveas(gcf,[filename '_LengthDistribution.fig'])
saveas(gcf,[filename '_LengthDistribution.png'],'png')

CorAngles(CorAngles(:,:)==0) = NaN;
CorCos(CorCos(:,:)==0) = NaN;
% figure; boxplot(CorCos)

Angles_Av = nanmean(CorAngles);
CosAv = nanmean(CorCos);
CosStd=nanstd(CorCos);
LogCos = log(CosAv);
X = 1:samplingdistance:size(CosAv,2)*samplingdistance;

%% Fit of log cos(theta)
[xData, yData] = prepareCurveData( X, LogCos );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedPoints = xData > P90;
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf 0];
opts.Upper = [0 0];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure;
h = plot( fitresult, xData, yData, excludedPoints, 'predobs' );
legend( h, 'LogCos vs. X', 'Excluded LogCos vs. X', 'linear Fit', 'Lower bounds (95% CI)', 'Upper bounds (95% CI)', 'Location', 'NorthEast' );
% Label axes
xlabel('Distance in nm')
ylabel('log < cos \theta >')
xlim([0 400])
ylim([-0.3 0.1])
grid on
hold off


%% 
Persistancelength = -1/fitresult.p1; % Persistence length P =  -s / ln (cos(theta));
% Persistancelength_std = -s/log(Cos_mean+Cos_std);
Rsquare=gof.rsquare;

disp(['Persistancelength (' num2str(s) ' nm spacing) = ' num2str(Persistancelength/1000) ' um   R^ = ' num2str(Rsquare) ])

data(1,5)= Persistancelength;
data(1,6)= Rsquare;
data(1,7)= P90; % points up to 90% of filament-lengths histogram included for fit 

saveas(gcf,[filename '_PLength.fig'])
saveas(gcf,[filename '_PLength.png'],'png')

csvwrite('mHtt_Data.csv',data);
% tom_emwrite('Htt_data.em',data);



% % Histogram cos(theta)
% % subplot(1,2,1,'Parent',figure); hold on
% a = axes('Parent',figure,'LineWidth',2,'FontName','Arial'); hold on
% h = histogram(CorCos,50);
% title('Distribution of cos(theta)');    
% xlabel('cos(theta)','fontsize',10,'fontweight','b'); 
% ylabel('Frequency','fontsize',10,'fontweight','b');
% set(a, 'FontSize', 10,'fontweight','b');
% set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 10, 'FontName', 'Arial');
%legend('cos(theta)','Location','NorthEast');


% % Histogram theta
% % subplot(1,2,2,'Parent',figure); hold on
% figure; a = axes;                                                    
% h = histogram(Angles_Av,50); 
% title('Distribution of theta');    
% xlabel('cos(theta)','fontsize',10,'fontweight','b'); 
% ylabel('Frequency','fontsize',10,'fontweight','b');
% set(a, 'FontSize', 10,'fontweight','b');
% set(gca, 'Box', 'on', 'LineWidth', 2, 'FontSize', 10, 'FontName', 'Arial');
% %legend('cos(theta)','Location','NorthEast');

%% test fitting each filament

LogCos = log(CorCos);
X = 1:samplingdistance:size(LogCos,2)*samplingdistance;
