%% This code analyzes minimal stimulation experiments (Figure 2)
% Results generated from this code are used to analyze Figure 2 data 

% Upload your data as many trials of synaptic stimulation failures or
% successes
% data = [your_minimal_stim_data];

%%  Create time vector 
Fs = 20000;                % Fs = 20 kHz = 20,000 points/ sec sampling rate 
t=((1:length(data))/Fs)';

% .........leak subtraction ..........
leak =abs(mean(data(Fs*.25:Fs*.65,:)));
noleak = data - repmat(leak,size(data,1),1); 
noleak=data+leak;

%% get input size, threshold of success vs failure from hist  

inputT = Fs*1.002:Fs*1.0157; % this is the time window over which the input occurs
% adjust it based on when you stimualte the cell 

inputs = abs(min(noleak(inputT,:))); 

% use a histogram to set success/failure threshold....
bw=100; %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
binvector = 0: bw :(max(inputs)+20);  % bin start: step: end
figure,histogram(inputs,binvector,'LineWidth',2.5,'FaceColor',...
    [0.8,0.8,0.8]);
box off
ylabel('Events/bin','FontName','Arial');
xlabel('IPSC size (pA)','FontName','Arial');


%% plot input size per trial 
% ****************** CHANGE THRESHOLD PER CELL ******************

%%%%%%%%%%%%%%%%%
threshold = 500; % this is based off the histogram generated above 
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% 
ipscT = Fs*1.0025; % this is the peak time of the IPSC 
%%%%%%%%%%%%%%%%%
numTrials=(1:1:size(inputs'));

trialColor=zeros(length(numTrials),3);
Ipeak= -noleak(ipscT,:); % current at peak IPSC time point (IPSC time)
for i = 1 : length(Ipeak) % conditional coloring of trial by input size 
	if Ipeak(i) > threshold
		trialColor(i, :) = [0,0,0]; % Black
	else
		trialColor(i, :) = [0.5,0.5,0.5]; % Gray 
	end
end
figure,scatter(numTrials,Ipeak,50* ones(length(Ipeak), 1), trialColor,...
    'LineWidth',2.0')
hold on
xlabel('Trial #')
ylabel('IPSC size (pA)')
xlim([0 max(numTrials)+1])
ylim([-40 max(inputs)+10])


%% plot exp fit, and successes and failures pt 1
%%% Fit decay with a single exponential: 
baseT=Fs*1.0415;
decay=nanmean(noleak(:,inputs>threshold),2);
decay=decay(ipscT:baseT);

f=fit(t(ipscT:baseT),decay,'exp1');
[m,gof]=fit(t(ipscT:baseT),decay,'exp1')
%figure,plot(f,t(ipscT:baseT),decay) % plots just the decay with the fit 


%% plot exp fit, and successes and failures pt 2 
Td=1000*1./400 % ***ENTER b VALUE generated from previous section (this is the tau decay)

% plot successes and failures, with decay 
figure,plot(t,noleak(:,:),'Color',[0.65 0.65 0.65])
hold on 
plot(t,nanmean(noleak(:,Ipeak>threshold),2),'black','LineWidth',2.5)
hold on 
plot(t,nanmean(noleak(:,Ipeak<threshold),2),'red','LineWidth',2.5)
hold on 
plot(f,t(ipscT:baseT),decay)
xlabel('time (s)')
ylabel('IPSC (pA)')
xlim([0.995 1.03]) 
ylim([min(-inputs)-100 150])
box off

%% Plot IPSC sizes as conductances 

% input conductances
% Iinput = g * Vhold 
Vhold=0.040;%%%%%%
%%%%%%%%%%%%%%%%%
g = (inputs*(1*10.^-12))./Vhold; % Units: pA / mV = S
g = (g)*(1*10.^9); % S --> nS

colI=transpose(sort(inputs));
colG=transpose(sort(g));
sizes=horzcat(colI,colG);

ipsc=inputs(inputs>threshold);
gipsc=ipsc*10.^-3./Vhold;
aSing=[mean(ipsc),std(ipsc),size(ipsc,2),std(ipsc)./sqrt(size(ipsc,2))...
    ,size(inputs,2),Td,mean(gipsc),std(gipsc),std(gipsc)./sqrt(size(gipsc,2))];

gBinVector = 0: bw./0.03*(1*10.^-3):(max(inputs)+20)./0.03*(1*10.^-3);
figure,histogram(g,gBinVector,'LineWidth',1.5,'FaceColor',[.8,.8,.8])
box off
ylabel('Counts/bin','FontName','Arial');
xlabel('Conductance (nS)','FontName','Arial');

