options =psoset('waitbar','off','maxfunevals',1000,'algorithmtype','synchronous','plot','off','Iterationsnoimprovement',10,'display','off','boundsmethod','bounce','modificationmethod','dynamic','numparticles',100,'inertiaReductionfraction',0.05,'velocityreductionfraction',0.05,'initializationmethod','rand','limitmaxvelocity',0)
options.TolX=0;
options.TolFun=0;
options.TolCen=0;
ub=100*ones(1,2); %ub=1000*ones(1,16);
lb=-ub;
runs = 50
wb = waitbar(0, 'Percentage Done');
for i=1:runs
    [result(i).X result(i).FVAL result(i).ExitFlag] = pso(@H1,[],lb,ub,options);
    RepeatNum = i;
    waitbar(i/runs,wb)
end

min = 0;
 
% for i=1:runs
%     if min >= result(i).FVAL,
%         min = result(i).FVAL;
%     end
% end

mincounter = 0;
for i=1:runs
    fprintf('Run: %d Value: %d\n',i,result(i).FVAL);
    fprintf('Position: %d\n', result(i).X);
    if min + 0.001 > result(i).FVAL,
        mincounter = mincounter + 1;
    end
end
min = min
percentage = mincounter/runs *100
close(wb)

% the following code will only work when the exitflag is returning the
% number of functions used in testing mode.
sum=0;
for i=1:runs,
    sum=sum+result(i).ExitFlag;
end
AvgFunEvalsUsed=sum/runs

