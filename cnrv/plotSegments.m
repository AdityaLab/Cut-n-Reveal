%plot timeseries with given segments
% you need to provide the output directory, original vis file and segment
% file
function [] = plotSegments(resultsoutputdir,visfile, segfile,clusNo)
paths2 = genpath('libs/linspecer');
addpath(paths2);
Y = csvread(visfile,1);
Y = Y';
segments = csvread(segfile);
figNameV = strcat('oscV_',mat2str(clusNo));
figurefilenameV=strcat(resultsoutputdir,fignameV,'.pdf');
figure;
N=250;
tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
%colorVar=['r','g','b','c','m','y','k'];
hold on
for j=1:size(Y,1)
    coloridx=mod(j,N);
    if coloridx==0
       coloridx=N;
    end
    plot(Y(j,:),'color',tempvar(coloridx,:),'Linewidth',1.3); %plot county timeseries.
end
xlim([0,size(Y,2)]); 
y1=get(gca,'ylim');
x1=get(gca,'xlim');
%clustercolor=['b','y','r','k','g'];
%clustercolor=['y','y','y','y','y'];
for j=1:size(segmentindicesV,2)
      line([segmentindicesV(j),segmentindicesV(j)],y1,'Color','k');
end    
set([gca],'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf,figurefilenameV);
print(figurefilenameV,'-dpdf');
hold off                                 
end