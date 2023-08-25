tic
data = [];
% x = [];
% data = gpuArray(x); while running onto GPU use this

channels=51100;
x=randi([1 200],1,channels);
y=randi([1 200],1,channels);
% x=100;
% y=100;

anntenas=2;

iterations = 3;
isFirstIteration=true;
parfor i = 1:channels
    csi_val=[];
    disp(" started simulaton for channel : %s",num2str(i));
    cfgHE = wlanHESUConfig;
    cfgHE.ChannelBandwidth = 'CBW20';  % Channel bandwidth
    cfgHE.NumSpaceTimeStreams = anntenas;     % Number of space-time streams
    cfgHE.NumTransmitAntennas = anntenas;     % Number of transmit antennas
    cfgHE.APEPLength = 1e3;            % Payload length in bytes
    cfgHE.ExtendedRange = false;       % Do not use extended range format
    cfgHE.Upper106ToneRU = false;      % Do not use upper 106 tone RU
    cfgHE.PreHESpatialMapping = false; % Spatial mapping of pre-HE fields
%     cfgHE.GuardInterval = 0.8;         % Guard interval duration
    cfgHE.HELTFType = 4;               % HE-LTF compression mode
    cfgHE.ChannelCoding = 'LDPC';      % Channel coding
%     cfgHE.MCS = 3;                     % Modulation and coding scheme
    
    fs = wlanSampleRate(cfgHE);
    chanBW = cfgHE.ChannelBandwidth;
    TimeSamples=50;
    cfgModel = winner2.wimparset;
    cfgModel.NumTimeSamples = TimeSamples;     % Time samples for CSI collection
    cfgModel.IntraClusterDsUsed = "no";   % No cluster splitting
%     cfgModel.SampleDensity = 200000;         % For lower sample rate
    cfgModel.PathLossModelUsed = "no";   % Turn off path loss
    cfgModel.ShadowingModelUsed = "no";  % Turn off shadowing
    cfgModel.SampleRate=fs;
%       cfgModel.SampleDensity = 20;
    cfgModel.RandomSeed = 10; % For repeatability
    % {3=B1-Urban micro-cell , 4=B2-Bad Urban micro-cell,10=C1-Suburban, 11=C2 - urban macro cell,14=D1-Rural macro-cell, 15=D2a}.
    scenario = 11;
    % finding distance between user and base station < 20 then skip that channel%
    distance = sqrt((x(i)-100)^2+(y(i)-100)^2);
    
    fprintf(' distance between AP and user is = %d',distance);
    if distance < 20
        continue;

    end
    
    cfgLayout = createUsersLayout(x(i),y(i),scenario,anntenas);
    
    winChannel = comm.WINNER2Channel(cfgModel,cfgLayout);
%     chanInfo = info(winChannel);
    
        [CSI1,pathDelays,finalCond] = winner2.wim(winChannel.ModelConfig,winChannel.LayoutConfig);
    
    [CSI2,~,finalCond] = winner2.wim(winChannel.ModelConfig,winChannel.LayoutConfig,finalCond);
    
    CSI = cellfun(@(x,y) cat(4,x,y),CSI1,CSI2,'UniformOutput',false);
%     [H2,~,finalCond] = winner2.wim(cfgModel,cfgLayout,finalCond);
    
    %
    CSI_Conv = cell2mat(CSI);
    CSI_Conv = abs(CSI_Conv);
    
%     %below code snppet gives PDP
%     a=size(CSI_Conv);
%     delaysPower=zeros(1,a(3));
%     for k = 1:a(3)
%         delaysPower(i)=mean(CSI_Conv(1,1,k,:),'all');
%         
%     end
%     hold on
%     impz(delaysPower');
%     hold on
%     % end PDP
    
    Y=zeros(1,24);
    maxNumErrors = 1;   % The maximum number of packet errors at an SNR point
    maxNumPackets = 100; % The maximum number of packets at an SNR point
    
    packetErrorRate = zeros(1,iterations);

    
    count =1;
     
    results = [];
    % snr = mcs_scemes';
     iter = linspace(1,iterations,iterations);
     iter=iter';
    tbl = table(iter);
    counter = 1;
    startMcs = 0;
    endMcs = 11;
    once = true;
    for mcs_ = startMcs:endMcs
        
         %changing MCS
        cfgHE.MCS=mcs_;
        for GI = [3.2,0.8]
            cfgHE.GuardInterval = GI;
        
    
         
        for isnr = 1:iterations % iterations
        % Set random substream index per iteration to ensure that each
        % iteration uses a repeatable set of random numbers
        stream = RandStream('combRecursive','Seed',99);
        stream.Substream = isnr;
        RandStream.setGlobalStream(stream);
        
        
        % Indices to extract fields from the PPDU
        ind = wlanFieldIndices(cfgHE);
    
        % Get occupied subcarrier indices and OFDM parameters
        ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE);
    
        % Account for noise energy in nulls so the SNR is defined per
        % active subcarrier
    %     packetSNR = snr(isnr)-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
        
        % Loop to simulate multiple packets
        numPacketErrors = 0;
        numPkt = 1; % Index of packet transmitted
        
        while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
            % Generate a packet with random PSDU
            psduLength = getPSDULength(cfgHE); % PSDU length in bytes
            txPSDU = randi([0 1],psduLength*8,1);
            
            tx = wlanWaveformGenerator(txPSDU,cfgHE);
        
            % Add trailing zeros to allow for channel delay
            txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];
        
            
        %         
    %         reset(winChannel);
            [rx,pathgains] = winChannel(txPad);
    
             
            rx_converted = cell2mat(rx);
            rx=rx_converted;%important line
    
            % Packet detect and determine coarse packet offset
            coarsePktOffset = wlanPacketDetect(rx,chanBW);
            if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
                numPacketErrors = numPacketErrors+1;
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
        
            lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
            coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
            rx = frequencyOffset(rx,fs,-coarseFreqOff);
               % Extract the non-HT fields and determine fine packet offset
            nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);
        
            % Determine final packet offset
            pktOffset = coarsePktOffset+finePktOffset;
        
            % If packet detected outwith the range of expected delays from
            % the channel modeling; packet error
            if pktOffset>50
                numPacketErrors = numPacketErrors+1;
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
        
            % Extract L-LTF and perform fine frequency offset correction
            rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
            fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
            rx = frequencyOffset(rx,fs,-fineFreqOff);
        
            % HE-LTF demodulation and channel estimation
            rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
            heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfgHE);
            [chanEst,pilotEst] = wlanHELTFChannelEstimate(heltfDemod,cfgHE);
            
            % Data demodulate
            rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
            demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE);
        
            % Pilot phase tracking
            demodSym = wlanHETrackPilotError(demodSym,chanEst,cfgHE,'HE-Data');
        
            % Estimate noise power in HE fields
            nVarEst = heNoiseEstimate(demodSym(ofdmInfo.PilotIndices,:,:),pilotEst,cfgHE);
        
            % Extract data subcarriers from demodulated symbols and channel
            % estimate
            demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);
            chanEstData = chanEst(ofdmInfo.DataIndices,:,:);
        
            % Equalization and STBC combining
            [eqDataSym,csi] = heEqualizeCombine(demodDataSym,chanEstData,nVarEst,cfgHE);
        
            % Recover data
            rxPSDU = wlanHEDataBitRecover(eqDataSym,nVarEst,csi,cfgHE,'LDPCDecodingMethod','norm-min-sum');
        
            % Determine if any bits are in error, i.e. a packet error
            packetError = ~isequal(txPSDU,rxPSDU);
            numPacketErrors = numPacketErrors+packetError;
        %         disp(txPSDU);
    %         disp(" packet errors "+numPacketErrors);
            numPkt = numPkt+1;
        end
        
        % Calculate packet error rate (PER) at SNR point
        packetErrorRate(isnr) = numPacketErrors/(numPkt-1);
    %     results(mcs_)=packetErrorRate;
       
        disp(['MCS ' num2str(cfgHE.MCS) ','...
              ' Guard Band ' num2str(cfgHE.GuardInterval)...
              ' Iterarions ' num2str(isnr) ...
              ' completed after ' num2str(numPkt-1) ' packets,'...
              ' PER:' num2str(packetErrorRate(isnr))]);
        
        end
        %% Plot Packet Error Rate vs SNR
    %     tbl = table(snr,Output1,Output2);
        MCS=packetErrorRate';
        tbl = addvars(tbl,MCS);
        % calculate average packet error rate above all snr
        % compare with threshold 10%
        %
        if mean(MCS,'all') < 0.01
    
            Y(counter)=1;
        %
        else 
            Y(counter)=0;
        end
    
        if counter==1
            
            csi_val=abs(chanEst);
        end
        counter=counter+1;
        release(winChannel);
        end
         
        
     end   
%      arr=[]; % used for plotting data purpose
%      for i=0:(endMcs-startMcs)*2
%         if i == 0
%             arr=[arr "MCS"];
%         
%         else
%             arr=[arr "MCS_"+i];
%         end
%         
%     end
%     figure;
     head(tbl,3)
%     %semilogy(tbl,"snr_",arr);
%     plot(tbl,"iter",arr);
%     grid on
%     legend
%     xlabel('Iterations');
%     ylabel('PER');
%     title(sprintf('PER for HE Channel %s, scenario:Bad urban macro-cell',cfgHE.ChannelBandwidth));
      
    disp(Y);
%     CSI_Conv = abs(CSI_Conv);
    data=[data;{csi_val,Y,finalCond}];
end

save('Dataset_50K_C2.mat','data')
toc
% % passing scenarios for sample datapoint
% sampleTestData = getSampleTest(4);
% trainedNet = trainNetwork(data,layers,options);
% YPred = classify(trainedNet,sampleTestData(1));
% YValidation = sampleTestData(2);
% accuracy = mean(YPred == YValidation);
% disp(accuracy);
        
   
 
