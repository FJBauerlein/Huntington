%% Script to analyze local ER mobility around IBs
% Part 3 of 3

% Felix JB Baeuerlein

% Baeuerlein et al. Cell 2017


%% reduce List of ER-Movies to selected ones
Result = Radial;
if exist('select_List')==1
for s=1:size(select_List,1)
    if select_List(s,1)==0
    Result.Average_ER(:,s) = NaN;
    Result.Average_GFP(:,s) = NaN;
    Result.Average_Var(:,s) = NaN;
    Result.Average_Var_norm(:,s) = NaN;
    Result.Average_Var_normPx(:,s) = NaN;
    end
end
    disp('-- List of ER_Movies reduced by select_List --')
else
    disp('-- All input ER_Movies taken for calculation --')
end

%% calculating population parameters
for s=1:size(Result.Average_ER,2)
    Result.Average_ER(:,s) = (Result.Average_ER(:,s) - nanmin(Result.Average_ER(1:20,s))); % removes background, measured as lowest Intensity inside the IB
    Result.Average_ER(:,s) = Result.Average_ER(:,s)./nanmax(Result.Average_ER(20:40,s)); % scales the ER Intensity to 1 between 0-2um close to the IB
    Result.Average_GFP(:,s) = Result.Average_GFP(:,s) - nanmin(Result.Average_GFP(:,s)); % scales GFP to 1 
    Result.Average_GFP(:,s) = Result.Average_GFP(:,s)./nanmax(Result.Average_GFP(:,s)); % scales GFP to 1 
    Result.Average_Var_norm(:,s) = Result.Average_Var_norm(:,s) - nanmin(Result.Average_Var_norm(1:20,s)); % scales Variance to 1
    Result.Average_Var_normPx(:,s) = Result.Average_Var_normPx(:,s) - nanmin(Result.Average_Var_normPx(1:20,s)); % scales Variance to 1
end

Result.group.Radius = Radial.Radius;
Result.group.ER_norm = nanmean(Result.Average_ER,2);
Result.group.GFP_norm = nanmean(Result.Average_GFP,2);
Result.group.Var_norm = nanmean(Result.Average_Var_norm,2);
Result.group.Var_normPx = nanmean(Result.Average_Var_normPx,2);

% 95 CI
Result.group.ER_norm_CI(:,1) = nanmean(Result.Average_ER,2)+1.96*nanstd(Result.Average_ER,1,2)./sqrt(size(Result.Average_ER,2)-sum(isnan(Result.Average_ER),2));
Result.group.ER_norm_CI(:,2) = nanmean(Result.Average_ER,2)-1.96*nanstd(Result.Average_ER,1,2)./sqrt(size(Result.Average_ER,2)-sum(isnan(Result.Average_ER),2));
Result.group.GFP_norm_CI(:,1) = nanmean(Result.Average_GFP,2)+1.96*nanstd(Result.Average_GFP,1,2)./sqrt(size(Result.Average_GFP,2)-sum(isnan(Result.Average_GFP),2));
Result.group.GFP_norm_CI(:,2) = nanmean(Result.Average_GFP,2)-1.96*nanstd(Result.Average_GFP,1,2)./sqrt(size(Result.Average_GFP,2)-sum(isnan(Result.Average_GFP),2));
Result.group.Var_norm_CI(:,1) = nanmean(Result.Average_Var_norm,2)+1.96*nanstd(Result.Average_Var_norm,1,2)./sqrt(size(Result.Average_Var_norm,2)-sum(isnan(Result.Average_Var_norm),2));
Result.group.Var_norm_CI(:,2) = nanmean(Result.Average_Var_norm,2)-1.96*nanstd(Result.Average_Var_norm,1,2)./sqrt(size(Result.Average_Var_norm,2)-sum(isnan(Result.Average_Var_norm),2));
Result.group.Var_normPx_CI(:,1) = nanmean(Result.Average_Var_normPx,2)+1.96*nanstd(Result.Average_Var_normPx,1,2)./sqrt(size(Result.Average_Var_normPx,2)-sum(isnan(Result.Average_Var_normPx),2));
Result.group.Var_normPx_CI(:,2) = nanmean(Result.Average_Var_normPx,2)-1.96*nanstd(Result.Average_Var_normPx,1,2)./sqrt(size(Result.Average_Var_normPx,2)-sum(isnan(Result.Average_Var_normPx),2));


%% plot the Radial Averages
figure;

subplot(3,2,1)
title('Radial averages of Htt97Q-GFP')
for sp1=1:size(Result.Average_GFP,2)
    hold on
    plot(Result.group.Radius, Result.Average_GFP(:,sp1),'Color',[0.8 0.8 0.8])
end
plot(Result.group.Radius, Result.group.GFP_norm, 'g', 'LineWidth',2.5);
plot(Result.group.Radius, Result.group.GFP_norm_CI, 'g--', 'LineWidth',1.5);
ylim([0 1.2])
xlabel('Distance from IB surface [\mum]');
ylabel('normalized Intensity');

subplot(3,2,3)
title('Radial averages of ERmCherry')
for sp2=1:size(Result.Average_ER,2)
    hold on
    plot(Result.group.Radius, Result.Average_ER(:,sp2),'Color',[0.8 0.8 0.8])
end
plot(Result.group.Radius, Result.group.ER_norm, 'r', 'LineWidth',2.5);
plot(Result.group.Radius, Result.group.ER_norm_CI, 'r--', 'LineWidth',1.5);
ylim([0 1.2])
xlabel('Distance from IB surface [\mum]');
ylabel('normalized Intensity');

subplot(3,2,5)
title('Intensity-normalized Variance')
for sp3=1:size(Result.Average_Var_normPx,2)
    hold on
    plot(Result.group.Radius, Result.Average_Var_normPx(:,sp3),'Color',[0.8 0.8 0.8])
end
plot(Result.group.Radius, Result.group.Var_normPx, 'b', 'LineWidth',2.5);
plot(Result.group.Radius, Result.group.Var_normPx_CI, 'b--', 'LineWidth',1.5);
ylim([0 2])
xlabel('Distance from IB surface [\mum]');
ylabel('Variance');

subplot(3,2,[2,4,6])
hold on
% plot(Result.group.Radius, Result.group.ER_norm, 'r', 'LineWidth',3);
% plot(Result.group.Radius, Result.group.GFP_norm, 'g', 'LineWidth',3);
% plot(Result.group.Radius, Result.group.Var_normPx./max(Result.group.Var_normPx), 'b', 'LineWidth',3); 
% plot(Result.group.Radius, Result.group.Var_normPx, 'b', 'LineWidth',3); 
% % plot(Result.group.Radius, Result.group.Var_norm, 'c:', 'LineWidth',3);
ylim([-0.1 1.3])
% 95 CI
% plot(Result.group.Radius, Result.group.GFP_norm_CI, 'g--', 'LineWidth',2);
shadedErrorBar(Result.group.Radius, Result.group.GFP_norm, Result.group.GFP_norm_CI(:,1)-nanmean(Result.Average_GFP,2), {'g', 'LineWidth',2.5},1)
% plot(Result.group.Radius, Result.group.ER_norm_CI, 'r--', 'LineWidth',2);
shadedErrorBar(Result.group.Radius, Result.group.ER_norm, Result.group.ER_norm_CI(:,1)-nanmean(Result.Average_ER,2), {'r', 'LineWidth',2.5},1)
% plot(Result.group.Radius, Result.group.Var_normPx_CI, 'b--', 'LineWidth',2);
shadedErrorBar(Result.group.Radius, Result.group.Var_normPx./max(Result.group.Var_normPx), Result.group.Var_normPx_CI(:,1)-nanmean(Result.Average_Var_normPx,2), {'b', 'LineWidth',2.5},1)

xlabel('Distance from IB surface [\mu m]');
ylabel('Variance');
hold off

% paper main Figure
figure
set(gcf, 'defaultAxesColorOrder', [0 0 0 ; 0 0 0])
hold on
shadedErrorBar(Result.group.Radius, Result.group.GFP_norm, Result.group.GFP_norm_CI(:,1)-nanmean(Result.Average_GFP,2), {'g', 'LineWidth',2.5},1)
shadedErrorBar(Result.group.Radius, Result.group.ER_norm, Result.group.ER_norm_CI(:,1)-nanmean(Result.Average_ER,2), {'r', 'LineWidth',2.5},1)
shadedErrorBar(Result.group.Radius, Result.group.Var_normPx./max(Result.group.Var_normPx), Result.group.Var_normPx_CI(:,1)-nanmean(Result.Average_Var_normPx,2), {'b', 'LineWidth',2.5},1)
yyaxis left
xlabel('Distance from IB surface [\mum]');
ylabel('Normalized Intensity');
ylim([-0.05 1.2])

yyaxis right
ylabel('Normalized variance');
ylim([-0.05 1.2])

% paper Suppl Figure
figure;

subplot(3,1,1)
%title('Radial averages of Htt97Q-GFP')
for sp1=1:size(Result.Average_GFP,2)
    hold on
    plot(Result.group.Radius, Result.Average_GFP(:,sp1),'Color',[0.8 0.8 0.8])
end
plot(Result.group.Radius, Result.group.GFP_norm, 'g', 'LineWidth',2.5);
plot(Result.group.Radius, Result.group.GFP_norm_CI, 'g--', 'LineWidth',1.5);
ylim([0 1.2])
%xlabel('Distance from IB surface [\mum]');
ylabel('Normalized intensity');

subplot(3,1,2)
%title('Radial averages of ERmCherry')
for sp2=1:size(Result.Average_ER,2)
    hold on
    plot(Result.group.Radius, Result.Average_ER(:,sp2),'Color',[0.8 0.8 0.8])
end
plot(Result.group.Radius, Result.group.ER_norm, 'r', 'LineWidth',2.5);
plot(Result.group.Radius, Result.group.ER_norm_CI, 'r--', 'LineWidth',1.5);
ylim([0 1.2])
%xlabel('Distance from IB surface [\mum]');
ylabel('Normalized intensities');

subplot(3,1,3)
%title('Intensity-normalized Variance')
for sp3=1:size(Result.Average_Var_normPx,2)
    hold on
    plot(Result.group.Radius, Result.Average_Var_normPx(:,sp3),'Color',[0.8 0.8 0.8])
end
plot(Result.group.Radius, Result.group.Var_normPx./max(Result.group.Var_normPx), 'b', 'LineWidth',2.5);
plot(Result.group.Radius, Result.group.Var_normPx_CI./max(Result.group.Var_normPx), 'b--', 'LineWidth',1.5);
ylim([0 2])
xlabel('Distance from IB surface [\mum]');
ylabel('Variance');


% saveas(gcf,['Analysis/IB' num2str(IB_count,'%03i') '_' fname_series '_RadialAverage_IB_' num2str(m) '.png'])
hold off