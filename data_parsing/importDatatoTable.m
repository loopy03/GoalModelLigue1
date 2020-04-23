function [T ] = importDatatoTable( years , league )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Filenames = {} ;
for i=1:length(league)
    for j=1:length(years)
        Filenames{end+1} = [pwd,'\csv files\',league{i},'_',num2str(years(j)),'.csv'] ;
    end
end
        

nleagues = length(Filenames) ;
res = {} ;

for n=1:nleagues
    res2{n} = importliguefile(Filenames{n});
    res(length(res)+1:length(res)+length(res2{n}),:) = res2{n} ;
end

k=1 ;
while k<=length(res)
    if isempty(res{k,1})
        res(k,:) = [] ;
    end
    k=k+1 ;
 end  
T = res2table(res) ;
   


end

