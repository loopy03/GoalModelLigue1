clear all
close all
clc

%% Chargement des données
load('Pred_L1_2020.mat')

%% Points actuels
actPtsH = [14 27 17 24 24 35 19 28 19 26 32 20 28 20 37 22 29 18 24 9]' ;
actPtsA = [9 12 20 10 6 14 21 28 15 14 8 17 13 7 31 19 21 12 14 4]' ;
actPts = actPtsH+actPtsA ;

%% Initialisation
allTeams = unique([probT.HomeTeam;probT.AwayTeam]) ;
nTeams = length(allTeams) ;
nGames = height(probT) ;

%% Nombre de points espéré (éspérance mathématique)
expPtsH = zeros(nTeams,1) ;
expPtsA = expPtsH ;

for ii=1:nTeams
    idx = strcmp(probT.HomeTeam,allTeams(ii)) ;
    prob = [probT.pH(idx) probT.pD(idx) probT.pA(idx)] ;
    expPtsH(ii) = sum(3*prob(:,1)+1*prob(:,2)) ;
    
    idx = strcmp(probT.AwayTeam,allTeams(ii)) ;
    prob = [probT.pH(idx) probT.pD(idx) probT.pA(idx)] ;
    expPtsA(ii) = sum(3*prob(:,3)+1*prob(:,2)) ;
end

expPts = expPtsH+expPtsA ;


%% Simulations de la fin de saison
nSim = 100000 ;
res = zeros(nGames,nSim) ;

% SImulation des matchs réstants
for ii=1:nGames
    prob = [probT.pH(ii) probT.pD(ii) probT.pA(ii)] ;
    res(ii,:) = randsrc(1,nSim,[1 2 3; prob]) ;
end

%% Calcul des points obtenus par chaque équipe
supPtsH = zeros(nTeams,nSim) ;
supPtsA = supPtsH ;

for ii=1:nTeams
    idx = strcmp(probT.HomeTeam,allTeams(ii)) ;
    resTmp = res(idx,:) ;
    supPtsH(ii,:) = 3*sum(resTmp==1) + 1*sum(resTmp==2) ;
    
    idx = strcmp(probT.AwayTeam,allTeams(ii)) ;
    resTmp = res(idx,:) ;
    supPtsA(ii,:) = 3*sum(resTmp==3) + 1*sum(resTmp==2) ; 
end
supPts = supPtsH+supPtsA ;

PtsH = actPtsH+supPtsH ;
PtsA = actPtsA+supPtsA ;
Pts = actPts+supPts ;

%% Calcul du classement final de chaque équipe
[~,I] = sort(Pts,'descend') ;
rank = zeros(nTeams,nSim) ;
for ii=1:nTeams
    rank(ii,:) = find(I==ii)'-nTeams*[0:nSim-1] ;
end

mRank = mean(rank,2) ;
[~,I] = sort(mRank) ;
   
%% figure histogramme classement
figure()
hold on
edges = 1:21 ;
k = 1 ;
brank = rank ; % Pour pouvoir avoir la barre la plus haute en couleur rouge

while k<=nTeams
    subplot(4,5,k)
    hold on
    
    % Histogramme tout rouge
    histogram(rank(I(k),:),'BinEdges',edges,'Normalization','probability',...
        'FaceColor',[0.8500 0.3250 0.0980]);
    N = histcounts(rank(I(k),:),edges) ;
     
    % Suppression de la barre la plus haute pour l'histogramme bleu
    [mMax,rMax] = max(N) ;
    ibrank = brank(I(k),:)==rMax ;
    brank(I(k),ibrank)=-100 ;
    
    % Histogramme tout bleu sans barre principale
    histogram(brank(I(k),:),'BinEdges',edges,'Normalization','probability',...
        'FaceColor',[0 0.4470 0.7410]);
    
    % écriture de la place la plus souvent atteinte
    at = rank(I(k),ibrank) ;
    text(at(1)+0.5,0.1+mMax/nSim,[num2str(at(1),'%4.0f')],'HorizontalAlignment','center','FontSize',8)
    
    
    % Mise en forme
    xlim([1 21])
    ylim([0 1])
    
    text(11,0.9,[allTeams{I(k)},' (',num2str(mRank(I(k)),'%4.1f'),')'],...
        'HorizontalAlignment','center','FontWeight','bold',...
        'FontSize',9) % Titre avec espérance entre parenthèses

    set(gca, ...
     'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XGrid'       , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.2:1, ...
    'XTick'       , [1 5 10 15 20], ...
    'LineWidth'   , 1         ); 
    ax=gca ;
    ax.XAxis.MinorTickValues = [1:20]+0.5;
    
    % Change the locations of the tick labels
    XTick = ax.XTick ;
    ax.XTick = XTick+0.5 ;
    % Change the tick labels themselves
    ax.XTickLabel = num2str(reshape(XTick,[],1),'%2.0f');     

    k = k+1 ;
end

annotation('textarrow',[0.5 0.5],[0.06 0.5],'string','Rang au classement', ...
      'HeadStyle','none','LineStyle', 'none', 'TextRotation',0,...
      'FontSize',12);
  
annotation('textarrow',[0.08 0.5],[0.6 0.6],'string','Probabilité', ...
      'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,...
      'FontSize',12);

%% figure histogramme points
figure()
hold on
k=1 ;
edges = 0:101 ;
bPts = Pts ;
while k<=nTeams
    subplot(4,5,k)
    hold on
    histogram(Pts(I(k),:),'BinEdges',edges,'Normalization','probability',...
        'FaceColor',[0.8500 0.3250 0.0980]);
    
    N = histcounts(Pts(I(k),:),edges) ; 
    [mMax,rMax] = max(N) ;
    ipts = Pts(I(k),:)==rMax-1 ;
    bPts(I(k),ipts)=800 ;
    
    histogram(bPts(I(k),:),'BinEdges',edges,'Normalization','probability',...
        'FaceColor',[0 0.4470 0.7410]);
    at = Pts(I(k),ipts) ;
    
    text(at(1)+0.5,0.025+mMax/nSim,[num2str(at(1),'%4.0f')],'HorizontalAlignment','center','FontSize',8)
    
    xlim([0 101])
    ylim([0 0.21])
    text(50,0.185,[allTeams{I(k)},' (',num2str(actPts(I(k))+expPts(I(k)),'%4.1f'),' pts)'],...
        'HorizontalAlignment','center','FontWeight','bold',...
        'FontSize',9)

    set(gca, ...
     'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'on'      , ...
    'XGrid'       , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.05:0.2, ...
    'XTick'       , [0:10:101], ...
    'LineWidth'   , 1         );  

    ax=gca ;
    ax.XAxis.MinorTickValues = [0:101]+0.5;
    
    % Change the locations of the tick labels
    XTick = ax.XTick ;
    ax.XTick = XTick+0.5 ;
    % Change the tick labels themselves
    ax.XTickLabel = num2str(reshape(XTick,[],1),'%3.0f');     
     
    k = k+1 ;
end
% annotation('textbox', [0.16, 0, 0.7, 0.1],'FitBoxToText','on',...
%     'string', 'Points au classement','LineStyle','none',...
%     'VerticalAlignment','middle','HorizontalAlignment','center')

annotation('textarrow',[0.5 0.5],[0.06 0.5],'string','Points au classement', ...
      'HeadStyle','none','LineStyle', 'none', 'TextRotation',0,...
      'FontSize',12);
  
annotation('textarrow',[0.08 0.5],[0.6 0.6],'string','Probabilité', ...
      'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,...
      'FontSize',12);


%% figures probabilités de classement
rankProbMat = zeros(nTeams) ;
for ii=1:nTeams
    for rr=1:nTeams
        rankProbMat(ii,rr) = length(find(rank(I(ii),:)==rr))/nSim ;
    end
end

figure()

imagesc(1-rankProbMat,'AlphaData',0.5); % Display correlation matrix as an image
set(gca, 'XTick', 1:nTeams); % center x-axis ticks on bins
set(gca, 'YTick', 1:nTeams); % center y-axis ticks on bins
set(gca, 'YTickLabel', allTeams(I)); % set y-axis labels
xlabel('Rang au classement')
set(gca, ...
    'TickDir'     , 'out'     , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );  
colormap('hot');

x = repmat(1:nTeams,nTeams,1); % generate x-coordinates
y = x'; % generate y-coordinates
% Generate Labels
for ii=1:nTeams
    for jj=1:nTeams
        t{ii,jj} = num2str(rankProbMat(ii,jj),'%4.2f') ;
    end
end
text(x(:), y(:), t, 'HorizontalAlignment', 'Center')

%% figure boxplot classement
figure()
boxplot([rank(I,:)]','Labels',allTeams(I),'LabelOrientation','inline',...
         'Notch','off','OutlierSize',4);
set(gca,'YDir','reverse')

ylim([1 20])
ylabel('Positions au classement')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'YGrid'       , 'on'      , ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 1:20,...
    'xtick',1:20,...
    'xticklabel',allTeams(I),...
    'LineWidth',1); 
 xtickangle(90)
 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',1); % Set width

%% figure boxplot points
figure()
hold on

boxplot([Pts(I,:)]','Labels',allTeams(I),'LabelOrientation','inline',...
         'Notch','off','OutlierSize',4,'Orientation','horizontal');
set(gca,'YDir','reverse')

xlim([0 100])
xlabel('Points au classement')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.0 .0] , ...
    'XGrid'       , 'on'      , ...
    'XMinorTick'  , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'XTick'       , 0:10:100,...
    'ytick',1:20,...
    'LineWidth',1); 
 
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',1); % Set width*

xax = get(gca,'XAxis'); % xax = get(ax,'XAxis'); 
set(xax,'TickLength',[.02 .02])

xline(12.5,'LineWidth',1)
xline(4,'LineWidth',1)
xline(100,'LineWidth',1)
 
med = median(Pts,2) ;
for k=1:nTeams
    text(4.5,k,allTeams(I(k)),'FontWeight','bold')
    text(3.5,k,num2str(k),'HorizontalAlignment', 'Right','FontWeight','bold')
    yline(k-0.5,'LineWidth',1)
    text(med(I(k))+0.2,k-0.04,num2str(med(I(k))),'FontSize',8,...
        'Position',[0 0.5 0])
end


%% Figures comparaison places en LDC (Marseille vs Lille vs Rennes)
figure()
hold on
edges = 50:90 ;
k = 1 ;
subplot(1,2,1)
hold on
while k<=nTeams
    
    if strcmp(allTeams{I(k)},'Marseille') || strcmp(allTeams{I(k)},'Lille') %|| strcmp(allTeams{I(k)},'Rennes')
    histogram(Pts(I(k),:),'BinEdges',edges,'Normalization','probability');
    end
    
    k = k+1 ;
end

    xlim([50 90])
    ylim([0 0.12])
    xlabel('Points au classement')
    ylabel('Probabilité')
    
    legend('Marseille','Lille')

    set(gca, ...
     'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XGrid'       , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.01:0.12, ...
    'XTick'       , [0:5:100], ...
    'LineWidth'   , 1         );  

    ax=gca ;
    ax.XAxis.MinorTickValues = [0:101]+0.5;
    
    % Change the locations of the tick labels
    XTick = ax.XTick ;
    ax.XTick = XTick+0.5 ;
    % Change the tick labels themselves
    ax.XTickLabel = num2str(reshape(XTick,[],1),'%3.0f');   
    
subplot(1,2,2)
hold on
k = 1 ;
idx=[] ;
while k<=nTeams
   
    if strcmp(allTeams{I(k)},'Marseille') || strcmp(allTeams{I(k)},'Lille') %|| strcmp(allTeams{I(k)},'Rennes')
        idx=[idx k] ;
    end
    k = k+1 ;
end

boxplot([Pts(I(idx),:)]','Labels',allTeams(I(idx)),'LabelOrientation','inline',...
         'Notch','off','OutlierSize',4,'Orientation','horizontal');
set(gca,'YDir','reverse')

xlim([50 90])
xlabel('Points au classement')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XGrid'       , 'on'      , ...
    'XMinorTick'  , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'XTick'       , 0:5:100,...
    'ytick',1:20,...
    'LineWidth',1); 

a = get(get(gca,'children'),'children') ;   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idxy=strcmpi(t,'box');  % Find Box objects
boxes=a(idxy);          % Get the children you need
set(boxes,'linewidth',1); % Set width*

    for k=1:length(idx)
    text(48,k,allTeams(I(idx(k))),'FontWeight','bold','HorizontalAlignment', 'Right')
    end

%% Figures comparaison places en LDC (Lille vs Rennes)
figure()
hold on
edges = 50:90 ;
k = 1 ;
subplot(1,2,1)
hold on
while k<=nTeams
    
    if strcmp(allTeams{I(k)},'Lille') || strcmp(allTeams{I(k)},'Rennes') %|| strcmp(allTeams{I(k)},'Rennes')
    histogram(Pts(I(k),:),'BinEdges',edges,'Normalization','probability');
    end
    
    k = k+1 ;
end

    xlim([50 80])
    ylim([0 0.12])
    xlabel('Points au classement')
    ylabel('Probabilité')
    
    legend('Lille','Rennes')
    
    %text(11,0.9,[allTeams{I(k)},' (',num2str(mRank(I(k)),'%4.1f'),')'],...
     %   'HorizontalAlignment','center','FontWeight','bold')
    %title([allTeams{I(k)},' (',num2str(mRank(I(k)),'%4.1f'),')'])

    set(gca, ...
     'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XGrid'       , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.01:0.12, ...
    'XTick'       , [0:5:100], ...
    'LineWidth'   , 1         );  

    ax=gca ;
    ax.XAxis.MinorTickValues = [0:101]+0.5;
    
    % Change the locations of the tick labels
    XTick = ax.XTick ;
    ax.XTick = XTick+0.5 ;
    % Change the tick labels themselves
    ax.XTickLabel = num2str(reshape(XTick,[],1),'%3.0f');   
    
subplot(1,2,2)
hold on
k = 1 ;
idx=[] ;
while k<=nTeams
   
    if strcmp(allTeams{I(k)},'Lille') || strcmp(allTeams{I(k)},'Rennes') %|| strcmp(allTeams{I(k)},'Rennes')
        idx=[idx k] ;
    end
    k = k+1 ;
end

boxplot([Pts(I(idx),:)]','Labels',allTeams(I(idx)),'LabelOrientation','inline',...
         'Notch','off','OutlierSize',4,'Orientation','horizontal');
set(gca,'YDir','reverse')

xlim([50 80])
xlabel('Points au classement')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XGrid'       , 'on'      , ...
    'XMinorTick'  , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'XTick'       , 0:5:100,...
    'ytick',1:20,...
    'LineWidth',1); 

a = get(get(gca,'children'),'children') ;   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idxy=strcmpi(t,'box');  % Find Box objects
boxes=a(idxy);          % Get the children you need
set(boxes,'linewidth',1); % Set width*

    for k=1:length(idx)
    text(48,k,allTeams(I(idx(k))),'FontWeight','bold','HorizontalAlignment', 'Right')
    end
    
%% Figures comparaison 18ème place (Nimes vs DIjon)
figure()
hold on
edges = 25:55 ;
k = 1 ;
subplot(1,2,1)
hold on
while k<=nTeams
    
    if strcmp(allTeams{I(k)},'Dijon') || strcmp(allTeams{I(k)},'Nimes') %|| strcmp(allTeams{I(k)},'Rennes')
    histogram(Pts(I(k),:),'BinEdges',edges,'Normalization','probability');
    end
    
    k = k+1 ;
end

    xlim([25 55])
    ylim([0 0.12])
    xlabel('Points au classement')
    ylabel('Probabilité')
    
    legend('Dijon','Nîmes')

    set(gca, ...
     'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XGrid'       , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.01:0.12, ...
    'XTick'       , [0:5:100], ...
    'LineWidth'   , 1         );  

    ax=gca ;
    ax.XAxis.MinorTickValues = [0:101]+0.5;
    
    % Change the locations of the tick labels
    XTick = ax.XTick ;
    ax.XTick = XTick+0.5 ;
    % Change the tick labels themselves
    ax.XTickLabel = num2str(reshape(XTick,[],1),'%3.0f');   
    
subplot(1,2,2)
hold on
k = 1 ;
%f=1;
idx=[] ;
while k<=nTeams
   
    if strcmp(allTeams{I(k)},'Dijon') || strcmp(allTeams{I(k)},'Nimes') %|| strcmp(allTeams{I(k)},'Rennes')
        idx=[idx k] ;
    end
    k = k+1 ;
end

boxplot([Pts(I(idx),:)]','Labels',allTeams(I(idx)),'LabelOrientation','inline',...
         'Notch','off','OutlierSize',4,'Orientation','horizontal');
set(gca,'YDir','reverse')

xlim([25 55])
xlabel('Points au classement')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XGrid'       , 'on'      , ...
    'XMinorTick'  , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'XTick'       , 0:5:100,...
    'ytick',1:20,...
    'LineWidth',1); 

a = get(get(gca,'children'),'children') ;   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idxy=strcmpi(t,'box');  % Find Box objects
boxes=a(idxy);          % Get the children you need
set(boxes,'linewidth',1); % Set width*

    for k=1:length(idx)
    text(23,k,allTeams(I(idx(k))),'FontWeight','bold','HorizontalAlignment', 'Right')
    end

    
        




            




    






       
    