%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_find_expl_exog3(resultsoutputdir,inputfile,segfile,laplacefile)
    paths2 = genpath('libs/linspecer');
    fprintf('results dir %s\n',resultsoutputdir);
    addpath(paths2);
    %resultsoutputdir = '../result/Harvey_exog/';
    %inputfile = '../data/Harvey_60min_sample_normalized.csv';
    %segfile = '../data/Harvey_60min_sample.csv';
    data = csvread(inputfile,1);
	%% Load Hurricane Data
    segments = csvread(segfile);
    c2 = size(data,2);
    cpV=[];
    sz = 30;%size(segments,2);
    for seg=1:sz
        r1 = 1;
        r2 = segments(seg);
        if seg>1
            r1=segments(seg-1);
        end
        X = data(r1:r2,1:c2);
        X=X';
        segV=find_expl_per_seg(resultsoutputdir,X,laplacefile,seg);
        if isempty(segV)==0
            cpV = [cpV, segV];
        end
        cpV = [cpV,segments(seg)];
    end
    seg = sz;
	r1 = segments(seg);
    r2 = segments(seg+1);%size(data,1);
    X = data(r1:r2,1:c2);
    X = X';
    segV=find_expl_per_seg(resultsoutputdir,X,laplacefile,seg+1);
    if isempty(segV)==0
        cpV = [cpV, segV];
    end
    cpV = [cpV,segments(seg+1)];
    cpV = sort(cpV);
    %write new cut points in file
    segmentfilenameV=strcat(resultsoutputdir,'segV','.csv'); % return this segment
    csvwrite(segmentfilenameV,cpV);
    %plot figure
    data = data';
    figurefilenameV=strcat(resultsoutputdir,'explV','.pdf');
    figure;
    N=10;
    tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
    %colorVar=['r','g','b','c','m','y','k'];
    sz2= 2500;%size(data,2)
    hold on
    for j=1:size(data,1)
        coloridx=mod(j,N);
        if coloridx==0
           coloridx=N;
        end
        plot(data(j,1:sz2),'color',tempvar(coloridx,:),'Linewidth',1.3); %plot county timeseries.
    end
    xlim([0,sz2]);            
    y1=get(gca,'ylim');
    x1=get(gca,'xlim');
    for j=1:size(cpV,2)
        contains = any(segments(:) == cpV(j));
        if contains==1
            line([cpV(j),cpV(j)],y1,'Color','r');
        else
            line([cpV(j),cpV(j)],y1,'Color','k');
        end
    end    
    set([gca],'FontSize', 18);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,figurefilenameV);
    print(figurefilenameV,'-dpdf');
    hold off
    close(gcf);
    fprintf('plotting figure\n');
    rmpath(paths2);    
end
