%% This code analyzes mini frequency and amplitude (Figure 1, Figure 1 - figure supplemnts 3 and 4)
% Results generated from this code are used to analyze the indicated figure data 

% Upload your data as many trials of mini recordings 
% data = [your_mini_data];
%%
clear eventWindow final Fs Hd i miniamps mult raster raw raw1 tempt times ui windo windowSize 

Fs = 20000;

% lowpass filter 
% N = filter order, Fc = cutoff freq for the point 6 dB below bandpass value 
d=fdesign.lowpass('N,Fc',50,200,Fs); % N,6dB,sampling rate
designmethods(d);
Hd = design(d);

raw1=data(4200:end,:); % exclude voltage step

s = size(data);
raw = filter(Hd,data);
     
raster.x = [];
raster.y = [];
raster.MF.x = [];
raster.MF.y = [];

%%
for i = 1:s(2) % # trials 
            
            
            temp = data(:,i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            window=0.0005; % DCN % set the size of the integration window 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            window = Fs*window;
            mult.a = 1:length(data);
            mult.a = repmat(mult.a',1,window);
            mult.P = mult.a + repmat(1:window,length(data),1); %integration of +window ms for each point
            mult.M = mult.a - 4*repmat(1:window,length(data),1);%integration of -window ms for each point
            mult.P(mult.P>length(data)) = length(data); % do not try to integrate after trace finishes
            mult.M(mult.M<1) = 1; % do not integrate prior to trace starting
            
            test = sum(temp(mult.P),2) - sum(temp(mult.M),2); % generation integration trace,  difference of integrated windows for each point
            events = double(diff(test(1:5:end)<-250)>0); % threshold integration, take every 5th point to avoid catching noise
            % Every 1st point counts too little, every 10th counts too much
            events = events';
            events(2:5,:) = 0;
            events = events(:); % resample
            
            MF = double(diff(test(1:5:end)>200)>0);
            MF = MF';
            MF(2:5,:) = 0;
            MF = MF(:);
            
            mult.MF = find(MF == 1);
            mult.MF = repmat(mult.MF',Fs*0.15,1);
            mult.MF = mult.MF + repmat((-0.05*Fs:1:(0.1*Fs-1))',1,size(mult.MF,2));
            mult.MF = mult.MF(:);
            mult.MF(mult.MF<1) = 1;
            mult.MF(mult.MF>length(MF)) = length(MF);
            
            MF = zeros(size(MF));
            MF(mult.MF) = 1;
            
            mult.events = find(events == 1); %get event times
            %sel = [0;diff(mult.events)>(Fs*0.0001)]; %remove events that are within 10 ms of each other
            %mult.events = mult.events(~~sel);
            events = zeros(size(events)); % rebuild event trace
            events(mult.events) = 1; % populate only with events spaced > 10 ms apart
            
            ui.IPSC(:,i) = events;
            
            x = 1:length(events);
            x = x(~~events);
            raster.x = [raster.x;x'];
            y = i*ones(length(x),1);
            raster.y = [raster.y;y];
            
            x = 1:length(MF);
            x = x(~~MF);
            raster.MF.x = [raster.MF.x;x'];
            y = i*ones(length(x),1);
            raster.MF.y = [raster.MF.y;y];
            
end
        
clear i s events mult test temp d x y MF
        
        %%
        
        windowSize = round(Fs*0.005);
        s = size(data);
        
        for i = 1:s(2)
            ui.bins{i} = [ui.IPSC(:,i);zeros(windowSize-rem(length(ui.IPSC(:,i)),windowSize),1)];
            ui.bins{i} = reshape(ui.bins{i},windowSize,length(ui.bins{i})/windowSize);
            ui.bins{i} = sum(ui.bins{i});
            ui.hist(:,i) = horzcat(ui.bins{i});
        end
        
        times.hist = (1:length(ui.hist(:,1)))*windowSize/Fs;
        
        clear s
        
        %% collect all events

        eventWindow = Fs*0.2;
        
        for i = 1:size(data,2);
            mult = find(ui.IPSC(:,i) == 1);
            mult = repmat(mult',eventWindow,1);
            mult = mult-eventWindow*0.2;
            mult = mult+repmat((1:eventWindow)',1,size(mult,2));
            mult(mult>length(data)) = length(data);
            mult(mult<1) = 1;
            temp = data(:,i);
            temp(1:Fs*0.3) = nan;
            temp(Fs:Fs*1.2) = nan;
            ui.events.sorted{i} = temp(mult);
        end
        
        ui.events.all = cell2mat(ui.events.sorted);
        sel = sum(isnan(ui.events.all),1);
        ui.events.all = ui.events.all(:,sel<100);
        ui.events.base = nanmean(ui.events.all(1:eventWindow*0.19,:),1);
        ui.events.rebound = nanmean(ui.events.all(eventWindow*0.5:eventWindow*0.7,:),1)-ui.events.base;
        ui.events.amp = min(ui.events.all((1:Fs*0.01)+eventWindow*0.2,:))-ui.events.base;
        
        times.events = (1:length(ui.events.all))/Fs - eventWindow*0.2/Fs;
        
        clear i sel

        %% compile into permanent summary matrix
        final.ui= ui;
        final.raster = raster;
        final.events = ui.events;
        

%%
% generate average event
for i = 1:size(final.events,2)
    ave.event(:,i) = nanmean(final.events(i).all,2);
    ave.base(i) = nanmean(ave.event(1:400,i),1);
    ave.event(:,i) = ave.event(:,i) - ave.base(i);
end
        
% calculate event frequency
for i = 1:size(final.events,2)
    ave.time(i) = numel(final.ui(i).IPSC);
    ave.noevents(i) = size(final.ui(i).events.amp,2);
    ave.freq(i) = ave.noevents(i)./ave.time(i)*Fs; 
end

miniamps=abs(ui.events.amp);
      
aveminifreq=ave.freq;
aveminiamp=mean(miniamps);

asumt=[aveminifreq,aveminiamp];
max(miniamps);

%% Check to make sure you're not counting noise 
%plot(ui.events.all(:,randi([size(ui.events.all:1,2)],1,1)))

events=(raster.x/Fs);

figure,%xline(1,'LineWidth',1.25)
hold on
% xline(1.01,'LineWidth',1.25)
plot(events,raster.y,'o','Color','k','LineWidth',1.5,'MarkerSize',2)
xlim([0.9 1.3])
xlim([0.95 1.2])

ax=gca;
ax.FontName='Arial';
ax.FontSize=12;
ax.LineWidth=1.5;
ax.TickDir='out';
ax.TickLength=[.025,.025];
ax.YTick = 0:size(data,2):size(data,2);
ax.XTick=.9:.1:1.3;
ylim([-5 size(data,2)+5])
box off
ylabel('Trial')
xlabel('Time (s)')


%% Verification
t=((1:length(data(:,1)))/Fs)'; %sec

xmin=1;
xmax=size(data,2);
randtrial=round(xmin+rand(1,1)*(xmax-xmin));
clear xmin xmax 

figure,plot(t,data(:,:)+0,'Color',[.5 .5 .5])
xlim([0.5 2])
%xlim([.95 1.3])
ylim([-400 0])
ylim([-1000 0])
hold on
plot(events(:,:),raster.y-101,'o','Color','k','LineWidth',1.5)

xlim([1.03 1.09])

%% Histogram and summary data 

%only include events 20 ms before and 30 ms after the stimulation if used 
peristim = events >0.8 & events <1.3 % burst 
%peristim = events >0.95 & events <1.04 % single stim
% ignore the above if you are not stimulating at all 

PSTevents=events(peristim)

binvector=0.8:0.001:1.3;
figure,histogram(PSTevents,binvector)
%'FaceAlpha',.1,'EdgeAlpha',.1)
hold on
% xline(1,'LineWidth',1.25)
% xline(1.01,'LineWidth',1.25)
% xline(1.02,'LineWidth',1.25)
% xline(1.03,'LineWidth',1.25)

box off
ylabel('Events/trial','FontName','Arial');
xlabel('time(ms)','FontName','Arial');
ax=gca;
histTick=0.95:.1:1.3;
ax.FontName='Arial';
ax.FontSize=14;
%ax.FontWeight='bold';
ax.LineWidth=1.5;
ax.TickDir='out';
ax.TickLength=[.025,.025];
ylim([0 70])




xlim([0.95 1.2])
xlim([0.95 1.04]) %single stim

[counts,bins]=histcounts(PSTevents,'BinWidth',.001) %
counts=counts';
bins=bins';



