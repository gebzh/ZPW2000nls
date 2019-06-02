%% Fig 5.4 the comparison of accuracies of different algs
load(strcat(pwd,'/data_files/Fig 5_4 simulation results.mat'));
load(strcat(pwd,'/data_files/Fig 5_4 simulation results patch.mat'));
f0_mse(:,4) = f0_mse_rRELAX; f1_mse(:,4) = f1_mse_rRELAX;
errorcount(:,4) = errorcount_rRELAX;

addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 6.3]);
% for the legend
semilogy(Tds, sqrt(f0_mse(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf0),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--k','MarkerSize',3.5);

% erasing
hold on;
semilogy(Tds, sqrt(f0_mse(:,1)),'-*w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-ow','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-^w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-dw','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf0),'-w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--w','MarkerSize',3.5);

% real plot
hold on; temp = sqrt(f0_mse(:,1))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'*k','MarkerSize',3.5);
hold on; temp = sqrt(f0_mse(:,2))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'ok','MarkerSize',3.5);
hold on; temp = sqrt(f0_mse(:,3))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'^k','MarkerSize',3.5);
hold on; temp = sqrt(f0_mse(:,4))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,1)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf0),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--k','MarkerSize',3.5);

axis([0.09 1 1e-4 max(sqrt(f0_mse(:,4)))]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.19 0.78 0.62];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontname{Times new roman}(a)\fontname{宋体}载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX',...
    '{\itCRB}_{\fontname{Times new roman}appf0}','{\itCRB}_{\fontname{Times new roman}simf0}'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 238 with the tag '% for chapter 5 plus' in columnlegend.m is uncommented 
box on;
figure; 
set(gcf,'unit','centimeters','position',[10 5 6 6.3]);
% for the legend
semilogy(Tds, sqrt(f1_mse(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf1),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--k','MarkerSize',3.5);

% erasing
hold on;
semilogy(Tds, sqrt(f1_mse(:,1)),'-*w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-ow','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-^w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-dw','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf1),'-w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--w','MarkerSize',3.5);

% real plot
hold on; temp = sqrt(f1_mse(:,1))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'*k','MarkerSize',3.5);
hold on; temp = sqrt(f1_mse(:,2))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'ok','MarkerSize',3.5);
hold on; temp = sqrt(f1_mse(:,3))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'^k','MarkerSize',3.5);
hold on; temp = sqrt(f1_mse(:,4))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,1)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf1),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--k','MarkerSize',3.5);

axis([0.09 1 1e-4 max(sqrt(f1_mse(:,2)))]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.19 0.78 0.62];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontname{Times new roman}(b)\fontname{宋体}低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX',...
    '{\itCRB}_{\fontname{Times new roman}appf0}','{\itCRB}_{\fontname{Times new roman}simf0}'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
box on;
% hint: note that the line 237 with the tag '% for chapter 5 plus' in columnlegend.m is uncommented 
%% Fig 5.5 the comparison of amounts of decoding errors of different algs
load(strcat(pwd,'/data_files/Fig 5_4 simulation results patch2.mat'));
load(strcat(pwd,'/data_files/Fig 5_4 simulation results.mat'));
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
temp1 = [Tds(1:12) Tds_periodo Tds(13:end)];
tmp = errorcount(:,1); temp2 = [tmp(1:12); errorcount_periodo(:,1);tmp(13:end)];
temp2 = temp2 + 1e-18;
semilogy(temp1, temp2,'-*k','MarkerSize',3.5);
hold on;
tmp = errorcount(:,2); temp2 = [tmp(1:12); errorcount_periodo(:,2);tmp(13:end)];
temp2 = temp2 + 1e-18;
semilogy(temp1, temp2,'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, errorcount(:,3)+1e-18,'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, errorcount(:,4)+1e-18,'-dk','MarkerSize',3.5);
axis([0.09 0.3 1e-1 max(max(errorcount))]);
axes = gca;
axes.XTick = [0.11 0.15 0.19 0.25 0.3];
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.19 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s'},'Fontsize',7.5);
ylabel('\fontname{宋体}译码错误次数','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
