%% This code analyzes train experiments (Figures 3 and 4)
% Results generated from this code are used to analyze Figure 3 and 4 data 

% import train data (stimulus artifact will be subtracted) 
% data = [your_train_data];
%% First, measure basic cell parameters (Cm, Ra, Rin, etc.) 

clear a Asumt ave ax cap1 EPSC1 EPSCs EPSCSS  Fs nostim obj PPR mult...
    raw rev Rmstep sel stim temps times tpulse traces Vstep done asumt 
Fs = 20000;

ave.trace = mean(data(:,1:1),2);
Fs = 20000;                % Fs = 20 kHz = 20,000 points/ sec
t=((1:length(data))/Fs)';

Vstep=0.0025; % 2.5 mV

%Mean leak per trial 
leaki=abs(mean(data(Fs*.05:Fs*.095,:))); 

all.leak =abs(mean(data(Fs*.25:Fs*.65,:)));
ave.leak=mean(all.leak);
%Max value of first cap. transient 
cap1=abs(min(data(Fs*0.098:Fs*0.105,:)));
Rs=0.003./(cap1-leaki)*10.^6;
mean(Rs)
%Rm
Rmstep=abs(mean(data(Fs*.11:Fs*.15,:)));
Rm=0.005./abs((Rmstep-leaki))*10.^6;

cap1=abs(min(data(Fs*0.048:Fs*0.065,:)));
all.Ra = -Vstep./(cap1-ave.leak)*10.^6;
Rmstep=abs(mean(data(Fs*.08:Fs*.17,:)));
all.Rin = .005./abs((Rmstep-ave.leak))*10.^6;
all.Cm =abs((sum(ave.trace(Fs*0.05:Fs*0.07)-mean(ave.trace(Fs*0.07:Fs*0.09)))/Fs)/Vstep); % pF

ave.Ra=mean(all.Ra)
ave.Rin=mean(all.Rin)
ave.Cm=mean(all.Cm)
%aaa=[ave.leak,ave.Rin,ave.Cm,ave.Ra,ave.Cm*ave.Rin./(10^6)*1000]




%%
all.EPSC = min(data(Fs*1.0015:Fs*1.005,:))-mean(data(Fs*0.95:Fs*0.99,:),1);

test = 0;

if test == 1;
    for i = 1:size(data,2)
    y = data(Fs*1.0015:Fs*1.01,:)-mean(data(Fs*0.95:Fs*0.99,:),1);
    [A,t] = min(y);
    y = y((t+3):end);
    x = (1:length(y))';
    taufit = fit(x,y,'exp1');
    temp = coeffvalues(taufit);
    all.tau1(i) = -1./temp(2)/Fs*1000;
    
    y = [data(Fs*1.504:Fs*1.512,i),data((Fs*1.514:Fs*1.522),i),data((Fs*1.524:Fs*1.532)+2,i),data((Fs*1.534:Fs*1.542)+3,i),data((Fs*1.544:Fs*1.552)+4,i)];
    y = nanmean(y,2)-mean(data(Fs*0.95:Fs*0.99,i),1);
    [A,t] = min(y);
    y = y((t+3):end);
    x = (1:length(y))';
    taufit = fit(x,y,'exp1');
    temp = coeffvalues(taufit);
    all.tau50(i) = -1./temp(2)/Fs*1000;
    end
end

%clear x y temp taufit A t  test i
%% remove stim artifact

Fs = 20000; % sampling rate
raw = data; % create copy of raw data to be processed
raw=raw(5000:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artthresh = 10;%%%%%%% threshold for where all stim artifacts pass through 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stim = diff(raw>artthresh)>0; % stims are those points where artthresh is exceeded
stim(1:Fs*0.5,:) = 0; % baseline value 
nostim = mode(sum(stim)); % number of stimuli 

%%% Only keep trials where # stimuli is what you expect
raw = raw(:,sum(stim) == nostim); 
stim = stim(:,sum(stim) == nostim);

mult = repmat(1:length(stim),size(stim,2),1)'; % matrix of stim times
mult = mult.*stim; % multiply by the stim times logical so all is 0 except stim times 
mult = mult(mult~=0) - 8; % collect stim times from stim time matrix

mult = reshape(mult,nostim,size(stim,2)); 
% reshape stim time matrix mult to Nostim rows, Notrials columns
% now mult is matrix of m = # stims, n = # trials and each value is the
% time of the stimulus

% reshape stim times into matrix where columns are traces
tpulse = (mult(2,:) - mult(2-1,:))/Fs; % get time of last pulse relative to end of train
% to get tpulse (=ISI), subtract time of 2nd stim to that of 1st stim, convert to ms

% add row numbers to each column to specify certain points in data matrix
mult = mult + repmat(0:size(raw,1):size(raw,1)*(size(raw,2)-1),size(mult,1),1);
% 0: size(raw,1) creates an index of all data point times for trial 1 
% size(raw,2)-1 is # trials - 1 
% size(raw,1)*(# trials - 1) gives total number of data points for all trials
% altogether, this gives each value of stim times in mult a unique index
% corresponding to the original data point 


% Create 3D matrices for artifact removal and EPSC selection 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artdur=0.00125; %%%%%%%%%%%% 0.00125 usec stim artifact duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

artifact = repmat(mult,[1 1 Fs*artdur]);
% 1st dimension encodes stim start time for each of the # of stimuli  
% 2nd dimension encodes trial / trace number 
% 3rd dimension encodes duration of stim time  

EPSClength = 0.0105; %0.0105
EPSCmult = repmat(mult,[1 1 (Fs*EPSClength+1)]); % why this +1? 

% Create stim time index for each trial and each stimulus 
admult = repmat((1:Fs*artdur) - 1,[size(artifact,2), 1, size(artifact,1)]);
% (1:Fs*artdur) = times for a single stim artifact starting from 0, not 1 
% size(artifact,2) = number of trials 
% size(artifact,1)] = number of stimuli 


admult = shiftdim(admult,2); % why do this? 

% artifact removal
raw(artifact+admult) = nan;

    
admult = repmat((-Fs*0.0005:Fs*0.01) - 1,[size(artifact,2), 1, size(artifact,1)]);

admult = shiftdim(admult,2);  

% EPSC selection
%sizead=size(admult);
%EPSCmult=EPSCmult(1:sizad(1),1:sizad(2),1:sizad(3));
EPSCs.all = raw(admult+EPSCmult);
% 1st dim = 211  
% 2nd dim is # stimuli 
% 3rd dim is # trials 
% so plot(EPSCs.all(:,:,1)) gives all EPSCs for all stimuli for trace 1 
EPSCs.all = shiftdim(EPSCs.all,2);

clear ans admult artifact EPSCmult mult ins rev EPSClength artdur artthresh

%% measure EPSC sizes

times.test = unique(round(tpulse*1000))/1000; % gives each unique stim frequency 
tpulse = round(tpulse*1000)/1000;
traces.sorted = [];
times.sorted = [];

% This will be EPSCs sorted by each unique stim frequency 
EPSCs.sorted = [];
for i = 1:length(times.test)
    EPSCs.sorted(:,:,i) = nanmean(EPSCs.all(:,:,tpulse == times.test(i)),3);
    test = cumsum(diff(isnan(EPSCs.sorted(:,1,i)))>0);
    if ~isempty(find(test >= 2))
        EPSCs.sorted(find(test >= 2),:,i) = nan;
    end
end

EPSCs.base = [shiftdim(nanmean(EPSCs.sorted(1:10,:,:)),1)];
EPSCs.size = shiftdim(min(EPSCs.sorted),1) - EPSCs.base; %%%%MEASURED FROM BASE


EPSCs.EPSC1 = mean(EPSCs.size(1,:)); % set EPSC1
EPSCs.ave = mean(EPSCs.size,2);
EPSCs.norm = EPSCs.ave(:)./EPSCs.ave(1);
% this is EPSCs.norm measured from start of stim artifact, not exp fit 



for i = 1:length(times.test)
    traces.sorted = [traces.sorted,raw(:,tpulse == times.test(i))];
    traces.ave(:,i) = mean(raw(:,tpulse == times.test(i)),2);
    traces.stimun(:,i) = mean(stim(:,tpulse == times.test(i)),2);
    traces.stimun(:,i) = traces.stimun(:,i) > 0.5;
    times.sorted = [times.sorted,tpulse(tpulse == times.test(i))];
    n(i) = sum(tpulse == times.test(i));
end

EPSCs.norm = EPSCs.size./repmat(EPSCs.size(1,:),size(EPSCs.size,1),1);
% this is each EPSC normalized to EPSC of the first EPSC of that frequency 


clear test
if exist('art','var')
    test.all = fieldnames(art);
    test.short = strfind(test.all,'short');
    test.long = strfind(test.all,'data');
    test.short = ~cellfun(@isempty,test.short);
    test.long = ~cellfun(@isempty,test.long);
else
    test.short = 0;
    test.long = 0;
end
% if artifact subtraction taken from single stim, subtract from all stims
% and reconstruct ave trace
if sum(test.short)>0
    % reconstruct ave trace with artifact subtraction
    traces.rec = traces.ave - nanmean(art.short(1:20));
   % for i = 1:size(traces.ave,2)
   %     temp =  EPSCs.sorted(:,:,i);
   %     traces.rec(:,i) = temp(:);
   % end
   
   for k = 1:size(traces.ave,2)
       for i = 1:sum(traces.stimun(:,k))
           ins = find(traces.stimun(:,k) == 1);
           %if length(traces.ave) < ins(end)+length(art.short)
           %    ins = ins(1:end-1);
           %end
           traces.rec(ins(i)-11:ins(i)+188,k) = traces.ave(ins(i)-11:ins(i)+188,k) - art.short(12:end);
       end
   end
end

%save(textdata{1})
traces.norm = traces.ave - repmat(mean(traces.ave(10000:12000,:),1),length(traces.ave),1);
%traces.norm = -traces.norm./repmat(min(traces.norm(Fs*1.5:Fs*1.505,:)),length(traces.norm),1);

if test.short == 1
traces.normrec = traces.rec - repmat(mean(traces.rec(10000:12000,:),1),length(traces.rec),1);
traces.normrec = -traces.norm./repmat(min(traces.norm(Fs*0.1:Fs*0.105,:)),length(traces.norm),1);
end

clear i test
%% measure EPSC sizes using single exp fit to IPSCn-1

if ~exist('art','var')
    for i = 1:length(times.test)
        EPSCs.sorted(:,:,i) = nanmean(EPSCs.all(:,:,tpulse == times.test(i)),3);
        test = cumsum(diff(isnan(EPSCs.sorted(:,1,i)))>0);
        if ~isempty(find(test >= 2))
            EPSCs.sorted(find(test >= 2),:,i) = nan;
        end
    end
end

 
 for j = 1:size(EPSCs.sorted,3)
     a = 1;
     
     for i = 1:size(EPSCs.sorted,2)
         [A t] = min(EPSCs.sorted(11:150,i,j));
         A = A - nanmean(EPSCs.sorted(1:11,1,j),1);
         t = t+10;
         x = (1:(length(EPSCs.sorted(:,i,j))-t))';
         y = EPSCs.sorted(t:length(EPSCs.sorted(:,i,j))-1,i,j) - nanmean(EPSCs.sorted(1:11,1,j),1);
         y = y(~isnan(y));
         x = x(~isnan(y));
         if length(y)>2
         taufit = fit(x,y,'exp1');
         e = coeffvalues(taufit);
         end
         %if j == 1
         %subplot(1,2,1),cla,plot(taufit,x,y),pause(1)
         %end
         
         if i == 1
             if abs((nanmean(EPSCs.sorted(1:11,1,j),1) - nanmean(EPSCs.sorted(1:11,i,j),1))) > 50
             EPSCs.expamp(1,j) = A;
             x = 1:1000;
             off = (e(1)*(exp(x.*e(2))))' + nanmean(EPSCs.sorted(1:11,1,j),1);
             off = off((length(y)):end);
             off = off(1:length(EPSCs.sorted(:,i+1,j)-1));
             %off = off - (mean(off(1:11)) - mean(EPSCs.sorted(1:11,i+1,j),1));
             test = EPSCs.sorted(:,i+1,j) - off;%+e(3)*(1-exp(x.*e(4))));
             EPSCs.expamp(i+1,j) = min(test(30:end));
             %if j == 1
             %    subplot(1,2,2),cla,plot(EPSCs.sorted(:,i+1,j)), hold on, plot(off);
             %   pause(1)
             %end
             else
                 EPSCs.expamp(1,j) = A;
                 EPSCs.expamp(i+1,j) = min(EPSCs.sorted(11:end,i+1,j)) - nanmean(EPSCs.sorted(1:11,i+1,j),1);
             end
         elseif i < size(EPSCs.sorted,2)
             if sum(times.test)<0
                 times.test = -times.test;
             end
             sel = times.test<0.04;
             if sel(j) > 0%abs((mean(EPSCs.sorted(1:11,1,j),1) - mean(EPSCs.sorted(1:11,i,j),1))) > 10
               x = 1:1000;
             off = (e(1)*(exp(x.*e(2))))' + nanmean(EPSCs.sorted(1:11,1,j),1);
             off = off((length(y)):end);
             off = off(1:length(EPSCs.sorted(:,i+1,j)-1));
             %off = off - (mean(off(1:11)) - mean(EPSCs.sorted(1:11,i+1,j),1));
             test = EPSCs.sorted(:,i+1,j) - off;%+e(3)*(1-exp(x.*e(4))));
             EPSCs.expamp(i+1,j) = min(test(30:end));
             else
                 EPSCs.expamp(i+1,j) = min(EPSCs.sorted(11:end,i+1,j)) - nanmean(EPSCs.sorted(1:11,i+1,j),1);
             end
             %if j == 1
            % subplot(1,2,2),cla,plot(EPSCs.sorted(:,i+1,j)), hold on, plot(off);
            % pause(1)
            %end
         end
         
     end

end      

EPSCs.normexp = EPSCs.expamp./repmat(EPSCs.expamp(1,:),size(EPSCs.expamp,1),1);
%this is EPSCs measured to extrapolated exp decay point 
%%  
rev = 1;
if rev == 1;
EPSCs.normexp = fliplr(EPSCs.normexp);
traces.ave = fliplr(traces.ave);
if isfield(traces,'rec')
traces.rec = fliplr(traces.rec);
end
%traces.rec = fliplr(traces.rec);
traces.norm = fliplr(traces.norm);
EPSCs.norm = fliplr(EPSCs.norm);
EPSCs.expamp = fliplr(EPSCs.expamp);
EPSCs.expamp = fliplr(EPSCs.expamp);
EPSCs.base = fliplr(EPSCs.base);
EPSCs.size = fliplr(EPSCs.size);
times.test = fliplr(times.test);
end
clear i j x y t taufit test e n off A addlater 


%%
EPSC1=EPSCs.ave(1);
PPR=EPSCs.ave(2)./EPSC1;
%EPSCSS=mean(EPSCs.norm(35:45)); % calculate steady-state synaptic strength
trials=EPSCs.norm;

figure
hold on
scatter(1:length(EPSCs.normexp),EPSCs.normexp,20,[0 0 0],'LineWidth',1.5)
%scatter(1:length(hthS),hthS,20,[0 0 0],'LineWidth',1.5)
%scatter(1:length(EPSCs.normexp),EPSCs.norm,20,[1 0 0],'LineWidth',1.5)
xlabel('Stimulus no.')
ylabel('IPSC (norm.)')
ylim([0 1.05])
box off
ax=gca;
ax.FontName='Arial';

ax.FontSize=14;
ax.LineWidth=1.5;
ax.TickDir='out';
ax.TickLength=[.025,.025];
%ax.XTick=0:50:150; %HTH_short
ax.XTick=0:50:200; % HTH
ax.YTick=0:.25:1;

%[m,gof]=fit((1:50)',EPSCs.norm(1:50),'exp1')

%%
xmin=1;
xmax=size(data,2);
randtrial=round(xmin+rand(1,1)*(xmax-xmin));
clear xmin xmax 

temps=(1:length(traces.norm))./Fs; 
figure,
plot(temps,traces.ave,'k','LineWidth',1.5) 
%plot(temps,traces.ave,'r','LineWidth',1.5) 
%plot(temps,traces.sorted(:,randtrial),'Color',[.5 .5 .5],'LineWidth',1.5)
box off
xlabel('time (s)')
ylabel('IPSC (pA)')

temps=temps';