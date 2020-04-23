clear all
close all
clc

%% Chargement de la librairie goalmodel
addpath([pwd,'/data_parsing'])
addpath([pwd,'/goalmodel_fit'])
addpath([pwd,'/goalmodel_pred'])

%% Chargement des données
years = 2018:2019 ;
league = {'F1','F2'} ;
T = importDatatoTable( years , league ) ;
load('csv files/F1_2019_upcoming.mat')

%% Initialisation 
% Date de départ de la simulation
initial_date = datetime('10/03/2020','InputFormat','dd/MM/yy') ;
initial_day = datenum(initial_date) ;

% Données connues
T_connue = T(T.Days<=initial_day,:) ;
team1_connue = T_connue.HomeTeam ;
team2_connue = T_connue.AwayTeam ;
goals1_connue = T_connue.HomeGoals ;
goals2_connue = T_connue.AwayGoals ;

% Données à prédire
T_inconnue = T_upcoming ;
team1_inconnue = cellstr(T_inconnue.HomeTeam) ;
team2_inconnue = cellstr(T_inconnue.AwayTeam) ;

% Propriétés du modèle
dc = true ; 
rs = false ;
model = 'poisson' ;
hfa = true ;
currentDate = datetime('21/04/2020','InputFormat','dd/MM/yy') ;
xi = 0.0065 ;
weights = weightsDC(T_connue.Dates,'xi',0.002,'currentDate',currentDate) ;

%% Goal Model computation
modelFit = goalModel(goals1_connue,goals2_connue,team1_connue,team2_connue,'dc',dc,'rs',rs,...
        'hfa',hfa,'model',model,'weights',weights) ;
    
%% predict Result 
prob = predictResult(modelFit,team1_inconnue,team2_inconnue) ;

%% predict ExpGoals
expG = predictExpg(modelFit,team1_inconnue,team2_inconnue) ;

%% Mise en forme
probT = table(team1_inconnue,team2_inconnue,prob(:,1),prob(:,2),prob(:,3),[expG.expg1]',[expG.expg2]') ;
probT.Properties.VariableNames = {'HomeTeam','AwayTeam','pH','pD','pA','expGH','expGA'} ;
