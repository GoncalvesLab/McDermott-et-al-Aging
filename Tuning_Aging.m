%This code is for 1x10 min long movies and
%with movement detection. It is based of Jake Jordan's version and ABC
%version

%subtract offset of position data (new baseline at zero)

clear all

filelistdata2 = dir('filelist') %location of first movie (any file will do)

MaxIterations=1; %set number of iterations

separateLaps=1; %set to (1) to analyze odd laps, (2) for even laps, (3) for first half of recording, (4) for last half of recording

TreadmillVersion=1; % Choose 0 for old treadmill (pre 01Jan2021) datasets A and B , choose 1 for new treadmill, dataset C

%set where each of the concatenated movies ends, so the transitions can be
%correctly evaluated

    levels = strfind(filelistdata2(1).folder, '\');

    mov1 = dir([filelistdata2(1).folder(1:levels(end)) 'CA1\Chan*_001*.tif']);
    %mov2 = dir([filelistdata2(1).folder(1:levels(end)) 'CA1_000\Chan*_001*.tif']);
    %mov3 = dir([filelistdata2(1).folder(1:levels(end)) 'CA1_001\Chan*_001*.tif']);
    
    movie1End=numel(mov1);
    %movie2End=numel(mov2);
    %movie3End=numel(mov3);
    movieEnds=movie1End;


%for files=1:size(filelistdata2,1)
    
    rpmdata=readNPY([filelistdata2(1).folder '\mset.npy']); %read RPM texture data
    rfiddata=readNPY([filelistdata2(1).folder '\tset.npy']); %read RFID texture data
    traces=readNPY([filelistdata2(1).folder '\cleaned_traces_data.npy']); %read deconvolved traces
    
    rfiddata=rfiddata'; %may be needed if loaded from .npy file
    rpmdata=rpmdata';   %may be needed if loaded from .npy file
    
    if TreadmillVersion==0
    
        detectionThreshold = 300;
    
    end

    if TreadmillVersion==1
    
        detectionThreshold = 700;
        rpmdata = rpmdata-1.372462080290498;
    end

    %First find transitions and resample RFID and RPM data to size of movies (traces)
    
    cd(filelistdata2(1).folder);
    
    %Find transitions in original 30kHz data as this is more robust
    
    transitions = findchangepts(rfiddata,'Statistic', 'mean', 'MinThreshold',detectionThreshold); %May need to tweak threshold, 300 works for 30kHz data
    
    %Plot transitions before trimming data
    
    hold off
    yvaltran=rfiddata(transitions);
    plot(rfiddata);
    hold on
    plot(transitions,yvaltran,'+r');
    savefig('transitions-pretrim_odd.fig');
    saveas(gcf,'transitions-pretrim_odd.jpg');
    
    %Then convert these points to timepoints in movie frame units
    
    movieTransitions=transitions.*(size(traces,2)/size(rfiddata,2));
    movieTransitions=ceil(movieTransitions);
    
    %Conversely, convert movieEnds to digitizer rpm/rfid time unids
    
    digitizerEnds = movieEnds./(size(traces,2)/size(rfiddata,2));
    digitizerEnds =round(digitizerEnds);
    if digitizerEnds(end)>size(rfiddata,2)
        digitizerEnds(end)=size(rfiddata,2);
    end

     %Remove initial part of each movie (except for first) until a zone transition is reached (as
    %cannot be sure of the position of the animal until that)
    
        
       %for movieBreak=(size(movieEnds,1)-1):-1:1
        %   transitionAfterNewMovie = movieTransitions-movieEnds(movieBreak);
        %   transitionAfterNewMovie(transitionAfterNewMovie < 0)=NaN; %set threshold as 0 frames but use ceil above because of rounding errors
        %   [~,transitionAfterNewMovie]=min(transitionAfterNewMovie);
        %   traces(:,movieEnds(movieBreak):movieTransitions(transitionAfterNewMovie))=[]; 
        %   rfiddata(digitizerEnds(movieBreak):transitions(transitionAfterNewMovie))=[];
        %   rpmdata(digitizerEnds(movieBreak):transitions(transitionAfterNewMovie))=[];       
       %end
       
     %Also remove first segment of first movie (until first transition is
     %reached)
     
     traces(:,1:movieTransitions(1))=[];
     rfiddata(1:transitions(1))=[];
     rpmdata(1:transitions(1))=[];

    %Now, make a new list of transitions in the trimmed traces (use 30KHz version) 
    
    transitions = findchangepts(rfiddata,'Statistic', 'mean', 'MinThreshold',detectionThreshold); %May need to tweak threshold, 300 works for 30kHz data
       
    
% Now begin Jake's code as we can assume that traces, rfiddata and rpmdata
% have all been trimmed to exlude areas that are not analyzeable.

rfidSession = zeros(1, size(rfiddata,2));
rfidSession(transitions) = 1;


%Find zones within RFID data and assign them a number


 % make a list of when the animal is in each zone...
    nzones = size(transitions,2)+1;
  
    zoneMeanVolt=zeros(nzones,1);
    beltZone=zeros(nzones,1);
    
    %...for this find mean voltage of each zone...
    
    for zone=1:1:nzones
        
        if (1<zone)&&(zone<nzones)
            zoneMeanVolt(zone)=mean(rfiddata(transitions(zone-1):transitions(zone)));
        elseif zone==1
            zoneMeanVolt(zone)=mean(rfiddata(1:transitions(zone)));
        elseif zone==nzones
            zoneMeanVolt(zone)=mean(rfiddata(transitions(zone-1):end));
        end    
    end
    
    %...and attribute each mean voltage to a texture
    
    if TreadmillVersion==0
        
        for zone=1:1:nzones
            if (round(zoneMeanVolt(zone)-0.60)==0)
                beltZone(zone)=1;
            elseif (round(zoneMeanVolt(zone)-1.8)==0)
                beltZone(zone)=2;
            elseif (round(zoneMeanVolt(zone)-2.3)==0)
                beltZone(zone)=3;
            elseif (round(zoneMeanVolt(zone)-1.2)==0)
                beltZone(zone)=4;

            end
        end
    end 
    
    if TreadmillVersion==1
        
          for zone=1:1:nzones
            if (round(zoneMeanVolt(zone)-2.4)==0)
                beltZone(zone)=1;
            elseif (round(zoneMeanVolt(zone)-0.60)==0)
                beltZone(zone)=2;
            elseif (round(zoneMeanVolt(zone)-1.2)==0)
                beltZone(zone)=3;
            elseif (round(zoneMeanVolt(zone)-1.8)==0)
                beltZone(zone)=4;

            end
          end
    end
    

    transitions(2,:) = transitions; % time points of transitions
    transitions(1,:) = beltZone(2:end);


% %cleaning raw encoder data and zeroing trace
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(rpmdata) %plot frameout
% [x,y] = ginput; %select 1 point to find noise threshold and hit return
% noiseThreshold = y(1,1); %first selected point
% rpmdata = rpmdata - noiseThreshold;
% rpmdata(rpmdata < 0) = 0; %remove everything below zero

d1 = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',0.15,'DesignMethod','butter');
rpmFiltered = filtfilt(d1,rpmdata);

mov = rpmFiltered;
%threshold rotation data to set baseline to zero
movSTDev = std(mov);
mov(mov < movSTDev*1.3) = 0; %this value was found just by looking at plot, seems far enough from baseline noise
mov = logical(mov);

%look for edges of running transients
edges=diff([0 mov]);
startruns = find(edges == 1);
endruns = find(edges == -1);

%Sometimes sessions start or end with a run in progress, so startruns and
%ednruns are not equal; need to remove the partial run
if size(startruns,2) ~= size(endruns,2)
    if startruns(1,1) < endruns (1,1)
        startruns = startruns(1:end-1); %removes a run session that ends after imaging ends 
    else
        endruns(1,1) = []; %removes a run session that begins prior to imaging session
    end
else
    startruns = startruns;
end

%find length of runs
lengthruns = endruns - startruns;

%define the start and end of stopped epochs (i.e. the start of a stop is
%the end of a run and vice-versa

startstops = endruns;
endstops= [startruns(2:end) 0];
lengthstops = endstops-startstops;

%find all stops lasting > 1sec (can adjust this time)

% startstops = startstops (lengthstops > frate*1);
% endstops = endstops (lengthstops > frate*1);

startstops = startstops (lengthstops > 15000);
endstops = endstops (lengthstops > 15000);

%create final trace as a flat y=1 line

movfinal=ones(size(mov));

%introduce "stopped" times as zeros in the trace

for i = 1:1:size(endstops,2)
movfinal(startstops(i):endstops(i))=0;
end

%now set beginning and end of trace (before first, after last run) to zero

movfinal(1:startruns(1))=0; 
movfinal(endruns(end):end)=0;

%Now we're left of running blips <1s (surrounded by stops > 1s)

%detect edges again to eliminate running blips <1s

edges = diff([0 movfinal]);
startruns = find(edges ==1);
endruns = find(edges ==-1);
lengthruns = endruns - startruns;

% startruns = startruns (lengthruns < frate*1);
% endruns = endruns (lengthruns < frate*1);

startruns = startruns (lengthruns < 15000);
endruns = endruns (lengthruns < 15000);

%set all "runs" shorter than 1s to zero

for i = 1:1:size(endruns,2)
movfinal(startruns(i):endruns(i))=0;
end

%% DO NOT DELETE
% Quality control for run selection performance
 movfinal(transitions(2,:))
 figure('units','normalized','outerposition',[0 0 1 1])
 subplot(2,1,1)
 plot(rpmdata,'r')
 hold on
 plot(movfinal)
 subplot(2,1,2)
 plot(rpmFiltered*5,'r')
 hold on
 plot(movfinal)
 for i = 1:1:size(transitions,2)
        if movfinal(transitions(2,i))==0
            disp("A transition point was cut off by the movement detection!")
        end
end
        


%%

%Resize treadmill data to trace files
movfinalResize = imresize(movfinal(1,:),[1,size(traces,2)]);
rfiddataResize = imresize(rfiddata(1,:),[1,size(traces,2)]);
rfidSessionResize = imresize(rfidSession(1,:),[1,size(traces,2)]);
rpmdataResize = imresize(rpmFiltered(1,:),[1,size(traces,2)]);

%Get address of all movement periods
movfinalResize = logical(movfinalResize);
rfidSessionResize(rfidSessionResize < 0.0001) = 0;
rfidSessionResize(rfidSessionResize > 0.0001) = 1;
movementIndex = find(movfinalResize == 1);

%Cut treadmill and traces data by movement index
movement = rpmdataResize(movementIndex);
rfiddataResize = rfiddataResize(movementIndex);
transitionsResize = rfidSessionResize(movementIndex);
tracesMove = traces(:,movementIndex); %resize trace files for movement periods

transitions(2,:) = find(diff(transitionsResize) == 1);

%Also remove last segment of  movie (after last transition is
%reached)
tracesMove(:,transitions(2,end):end)=[];
rfiddataResize(transitions(2,end):end)=[];
movement(transitions(2,end):end)=[];

%Also remove first segment of movie (until first transition is
%reached)
tracesMove(:,1:transitions(2,1))=[];
rfiddataResize(1:transitions(2,1))=[];
movement(1:transitions(2,1))=[];

transitions(2,:)= transitions(2,:) - (transitions(2,1)-1);
transitions(2,end) = size(movement,2) + 1;


% transitions(:,end) = [];
% transitions(:,1) = [];

%%
% firstTransition = transitions(2,1) - 1;
% lastTransition = transitions(2,end) + 1;
% 
% rfiddataResize(lastTransition:end) = [];
% movement(lastTransition:end) = [];
% tracesMove(lastTransition:end) = [];
% rfiddataResize(1:firstTransition) = [];
% movement(1:firstTransition) = [];
% tracesMove(1:firstTransition) = [];
% transitions(2,:) = transitions(2,:) - (transitions(2,1) - 1);

%Then make a list of when the animal is in each zone...


  
% zoneMeanVolt = zeros(nzones,1);

    
%...for this find mean voltage of each zone...

%     for zone = 1:1:nzones
%         
%         if (1 < zone) && (zone < nzones) %mean voltage of zones between the second and second to last transition
%             zoneMeanVolt(zone) = mean(rfiddataResize(transitionsFr(zone - 1):transitionsFr(zone)));
%         elseif zone == 1 %mean voltage of zone between first and second transition
%             zoneMeanVolt(zone) = mean(rfiddataResize(transitionsFr(zone)));
%         elseif zone == nzones %mean voltage of last transition
%             zoneMeanVolt(zone) = mean(rfiddataResize(transitionsFr(zone - 1):end));
%         end    
%     end
        
    
%     for zone = 1:1:nzones
%         zoneMeanVolt(zone) = mean(rfiddataResize(transitionsFr(zone):transitionsFr(zone + 1)));
%     end    
    
    %...and attribute each mean voltage to a texture
    
%     for zone = 1:1:nzones
%         
%         zoneVoltSieve = zeros(1,4);
%         zoneVoltSieve(1) = zoneMeanVolt(zone) - 1.20;
%         zoneVoltSieve(2) = zoneMeanVolt(zone) - 1.80;
%         zoneVoltSieve(3) = zoneMeanVolt(zone) - 2.30;
%         zoneVoltSieve(4) = zoneMeanVolt(zone) - 0.60;
%         
%         [~,beltZone(zone)] = min(abs(zoneVoltSieve));
%         
%     end
    
    %Create variable with start and end of each zone
    
    transitionsFr = transitions(2,:); %timestamps of all zone tag swipes; resized to traces and selected for movement periods
    beltZone = transitions(1,:)'; %IDs of all zone tag swipes; resized to traces and selected for movement periods
    beltZone(end) = []; %prevents code from trying to analyze past the last zone tage swipe
    nzones = size(transitionsFr,2) - 1; %number of texture zones being analyzed
    
    zoneStartEnd = zeros(nzones,2);
    
%     for zone = 1:1:nzones
%         
%         if (1 < zone) && (zone < nzones)
%             zoneStartEnd(zone,:) = [transitionsFr(zone - 1),(transitionsFr(zone) - 1)]; % subtract 1 so that border point isn't accounted twice
%         elseif zone == 1
%             zoneStartEnd(zone,:) = [1,(transitionsFr(zone) - 1)];
%         elseif zone == nzones
%             zoneStartEnd(zone,:) = [transitionsFr(zone - 1),size(rfiddataResize,2)];
%         end    
%     end

    for zone = 1:1:nzones
         zoneStartEnd(zone,1) = [transitionsFr(zone)];
         zoneEnd = zone +1;
         zoneStartEnd(zone,2) = [(transitionsFr(zoneEnd) - 1)]; % subtract 1 so that border point isn't accounted twice
    end
    
%     zoneStartEnd(1,:) = [];
%     zoneStartEnd(end,:) = [];
%     zoneStartEnd(end,end) = zoneStartEnd(end,end) - 1;
    %nzones = nzones -2;
    
    %Find distance for each trip through each texture
    for zone = 1:nzones
        distmeasure(zone) = sum(movement(zoneStartEnd(zone,1):zoneStartEnd(zone,2))); %Find distance for each trip through each texture
    end

%     beltZone(distmeasure < 2.2) = []; %delete partial zones
%     zoneStartEnd(distmeasure < 2.2,:) = [];
%     distmeasure(distmeasure < 2.2) = []; %delete partial zones
%     nzones = size(beltZone,1);
    
    %Count number of times that mouse travels through each of 4 textures

    textureone = zeros(sum(beltZone == 1),1);
    texturetwo = zeros(sum(beltZone == 2),1);
    texturethree = zeros(sum(beltZone == 3),1);
    texturefour = zeros(sum(beltZone == 4),1); 
    
    textureone = distmeasure(beltZone == 1); %distance traveled in texture one
    texturetwo = distmeasure(beltZone == 2); %distance traveled in texture two
    texturethree = distmeasure(beltZone == 3); %distance traveled in texture three
    texturefour = distmeasure(beltZone == 4);%distance traveled in texture four
    
    %Display zone errors for quality control
    percentError = [(max(textureone)-mean(textureone))/mean(textureone), (max(texturetwo)-mean(texturetwo))/mean(texturetwo), (max(texturethree)-mean(texturethree))/mean(texturethree), (max(texturefour)-mean(texturefour))/mean(texturefour)]
    meanZoneError = mean(percentError)
    
    
    %Find number of data points for each zone
    
    zonesize = zeros(nzones,1); %initialize variable
    
    for zone = 1:nzones
        zonesize(zone,1) = zoneStartEnd(zone,2) - zoneStartEnd(zone,1); % number of data points for each trip through each texture
    end
    
    zonesizeone = zonesize(beltZone == 1); %number of data points traveled in texture one
    zonesizetwo = zonesize(beltZone == 2); %number of data points traveled in texture two
    zonesizethree = zonesize(beltZone == 3); %number of data points traveled in texture three
    zonesizefour = zonesize(beltZone == 4);%number of data points traveled in texture four
    
    %Find cumulative distance along each zone, this needs to be a cell array
    %as each zone is a different size
    
    distancex{nzones,1} = []; %initialize cell array
    
    %Create x (distance) for each zone crossing
    
    for zone = 1:nzones
        distancex{zone} = abs(cumsum(movement(zoneStartEnd(zone,1):zoneStartEnd(zone,2)))); %integrate distance for each trip through each texture
    end
    
    %and organize this data by zones and laps
    
    xone = cell(size(textureone));
    newzone = 1;
    for zone = 1:1:nzones
        if beltZone(zone) == 1
            xone{newzone} = distancex{zone}./max(distancex{zone}); % x positions of data points in texture one
            xone{newzone} = xone{newzone}.*25;
            newzone = newzone+1;
        end
    end
    
    xtwo = cell(size(texturetwo));
    newzone = 1;
    for zone = 1:1:nzones
        if beltZone(zone) == 2
            xtwo{newzone} = distancex{zone}./max(distancex{zone}); % x positions of data points in texture two
            xtwo{newzone} = xtwo{newzone}.*25+25;
            newzone = newzone+1;
        end
    end
    
    xthree = cell(size(texturethree));
    newzone = 1;
    for zone = 1:1:nzones
        if beltZone(zone) == 3
            xthree{newzone} = distancex{zone}./max(distancex{zone}); % x positions of data points in texture three
            xthree{newzone} = xthree{newzone}.*25+50;
            newzone = newzone+1;
        end
    end
    
    xfour = cell(size(texturefour));
    newzone = 1;
    for zone = 1:1:nzones
        if beltZone(zone) == 4
            xfour{newzone} = distancex{zone}./max(distancex{zone}); % x positions of data points in texture four
            xfour{newzone} = xfour{newzone}.*25+75;
            newzone = newzone+1;
        end
    end

    %Initialize position variable
    AnimalPosition = 1;
    Lap = 1;
    
    %Initialize counters for each zone
    
    z1counter = 1;
    z2counter = 1;
    z3counter = 1;
    z4counter = 1;
    
    for zone = 1:1:nzones
        if beltZone(zone) == 1
            AnimalPosition = [AnimalPosition, xone{z1counter}];
            Lap = [Lap, ones(size(xone{z1counter}))*z1counter];
            z1counter = z1counter+1;
        elseif beltZone(zone) == 2
            AnimalPosition = [AnimalPosition, xtwo{z2counter}];
            Lap = [Lap, ones(size(xtwo{z2counter}))*z2counter];
            z2counter = z2counter+1;
        elseif beltZone(zone) == 3
            AnimalPosition = [AnimalPosition, xthree{z3counter}];
            Lap = [Lap, ones(size(xthree{z3counter}))*z3counter];
            z3counter = z3counter+1;
        elseif beltZone(zone) == 4
            AnimalPosition = [AnimalPosition, xfour{z4counter}];
            Lap = [Lap, ones(size(xfour{z4counter}))*z4counter];
            z4counter = z4counter+1;
        end
    end
    
    AnimalPosition(1) = [];
    Lap(1) = [];
    
    figure;
    plot(AnimalPosition);
    savefig('AnimalPosition_odd.fig');
%     saveas(gcf,'AnimalPosition.jpg');
            
    save('AnimalPosition_odd.mat', 'AnimalPosition');
%     save('Lap.mat', 'Lap');

    traceLog = logical(tracesMove); %binarize activity traces

    % Separate laps ? If different from 0, then only a subset of data will
    % be analyzed

    if separateLaps ~= 0

        [~,lapStart] = findpeaks(AnimalPosition,'MinPeakHeight', 99.9);
        lapStart=flip(lapStart);
    end
            
    if separateLaps == 1 %select odd laps only

        for point=1:2:size(lapStart,2)
            if (point+1) <= size(lapStart,2)
                AnimalPosition(lapStart(point+1):lapStart(point))=[];
            else
                 AnimalPosition(lapStart(point):1)=[];
            end
            
            if (point+1) <= size(lapStart,2)
               traceLog(:,lapStart(point+1):lapStart(point))=[];
            else
               traceLog(:,lapStart(point):1)=[];
            end       
        end
    end

    if separateLaps == 2 %select even laps only

        for point=2:2:size(lapStart,2)
            if (point+1) <= size(lapStart,2)
                AnimalPosition(lapStart(point+1):lapStart(point))=[];
            else
                 AnimalPosition(lapStart(point):1)=[];
            end
            
            if (point+1) <= size(lapStart,2)
               traceLog(:,lapStart(point+1):lapStart(point))=[];
            else
               traceLog(:,lapStart(point):1)=[];
            end      
        end
    end

    if separateLaps == 3 %select first half of trace only

        traceMidpoint=round(size(AnimalPosition,2)/2);
        AnimalPosition(traceMidpoint:end)=[];
        traceLog(:,traceMidpoint:end)=[];
    end
            
    if separateLaps == 4 %select second half of trace only

        traceMidpoint=round(size(AnimalPosition,2)/2);
        AnimalPosition(1:traceMidpoint)=[];
        traceLog(:,1:traceMidpoint)=[];
    end
   
    %Find time spent at each position
    
    edges = (0:1:100);
    polangle = linspace(0,2*pi,100); %use 100 bins for whole lap
    
    histTime = histcounts(AnimalPosition,edges);
    
    %Now find firing position for all neurons and build histogram
    
    
    
    
    %initialize output structure

    CellTuning{size(traceLog,1),8}=[];
    
    MaxIterations = 10000;
    
    for neuron = 1:1:size(traceLog,1)
        
        FirePos = AnimalPosition(traceLog(neuron,:));
        histLap = histcounts(FirePos,edges);
        histLap = histLap./histTime; %normalize to time spent at each position
        histLap = histLap/sum(histLap); % then normalize so sum = 1
        
        tuningphasor = 0;
        for i = 1:size(histLap,2)
            
            tuningphasor = tuningphasor + (histLap(i)*exp(1i*polangle(i)));
            
        end
        
        tuningidx = abs(tuningphasor);
        
        %now find if this is significant by bootstraping to create null
        %distribution
        
        randTuningDist = tuningShuffleNP(traceLog(neuron,:), MaxIterations, edges, histTime, AnimalPosition, polangle);
         
        SigThreshold = prctile(randTuningDist,95);
        SortedItrTuning = sort(randTuningDist);
        IndicesHigher = find(SortedItrTuning > tuningidx);
        if numel(IndicesHigher) ~= 0
    
            PVal = (MaxIterations-IndicesHigher(1))/MaxIterations;
        else
            PVal = 1/MaxIterations;
        end
        
        %Save information of each cell: the firing histogram (rho),the
        %angles (theta), the tuning index and the alpha value for 95th percentile
        
        CellTuning{neuron,1} = histLap;
        CellTuning{neuron,2} = polangle;
        CellTuning{neuron,3} = tuningidx;
        CellTuning{neuron,4} = SigThreshold; %threshold for alpha=0.05
        CellTuning{neuron,5} = [25, 50, 75, 100]; %borders between zones (angles)
        CellTuning{neuron,6} = PVal; %p-value of tuningidx
        CellTuning{neuron,7} = randTuningDist; %save null distribution
        CellTuning{neuron,8} = AnimalPosition; %save null distribution
        
    end
    
    save('CellTuning_odd.mat', 'CellTuning');
    
    
    %clearvars -except filelistmovtex filelistdata2 files MaxIterations %clears all vars except those needed to loop through different mice
    
    
    
    
%     clearvars -except beltData1 beltData2 imagingSessionXML cleanedTracesDir files MaxIterations CellTuning%clears all vars except those needed to loop through different mice
    

session1 = CellTuning; %Load CellTuning file as session1

normalizedMaps = [];
for n = 1:size(session1,1)
   peak_rate = max(session1{n,1}); %should change so that rasters are organized by center of mass; right now, we're organizing rasters by their peak average firing rate
   scaled = (1/peak_rate)*session1{n,1};
   normalizedMaps = [normalizedMaps; scaled]; 
end

cellPVals = [];
for n = 1:size(session1,1)
    pValue = session1{n,6};
    cellPVals = [cellPVals; pValue];
end

tunedCells = find(cellPVals < 0.05);
nonTunedCells = find(cellPVals >= 0.05);
tunedMap = normalizedMaps(tunedCells,:);
nonTunedMap = normalizedMaps(nonTunedCells,:);
[tunedMax, maxLoc] = max(tunedMap');
maxLoc = maxLoc';
maxLoc(1:size(maxLoc),2) = 1:size(maxLoc);
[~,idx] = sort(maxLoc(:,1)); % sort just the first column
sortLoc = maxLoc(idx,:);   % sort the whole matrix using the sort indices

tunedMapSorted =[];
for n = 1:size(tunedCells)
    map = tunedMap(sortLoc(n,2),:);
    tunedMapSorted = [tunedMapSorted; map];
end

allCellsMap = [tunedMapSorted; nonTunedMap];
figure;
heatmap(allCellsMap);
grid off
colormap(parula)
savefig('Raster_odd.fig');

    msgbox("Done")