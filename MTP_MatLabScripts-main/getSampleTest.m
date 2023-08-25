function f = getSampleTest(scenario)
    
    cfgHE = wlanHESUConfig;
    cfgHE.ChannelBandwidth = 'CBW20';  % Channel bandwidth
    cfgHE.NumSpaceTimeStreams = 2;     % Number of space-time streams
    cfgHE.NumTransmitAntennas = 2;     % Number of transmit antennas
    cfgHE.APEPLength = 1e3;            % Payload length in bytes
    cfgHE.ExtendedRange = false;       % Do not use extended range format
    cfgHE.Upper106ToneRU = false;      % Do not use upper 106 tone RU
    cfgHE.PreHESpatialMapping = false; % Spatial mapping of pre-HE fields
    cfgHE.GuardInterval = 0.8;         % Guard interval duration
    cfgHE.HELTFType = 4;               % HE-LTF compression mode
    cfgHE.ChannelCoding = 'LDPC';      % Channel coding
    cfgHE.MCS = 3;                     % Modulation and coding scheme
    
    fs = wlanSampleRate(cfgHE);
    chanBW = cfgHE.ChannelBandwidth;
    
    channels=1;
    x=randi([1 200],1,channels);
    y=randi([1 200],1,channels);
    TimeSamples = 1600;
    
    
    cfgModel = winner2.wimparset;
    cfgModel.NumTimeSamples = TimeSamples;     % Frame length
    cfgModel.IntraClusterDsUsed = "no";   % No cluster splitting
    cfgModel.SampleDensity = 2e5;         % For lower sample rate
    cfgModel.PathLossModelUsed = "yes";   % Turn on path loss
    cfgModel.ShadowingModelUsed = "yes";  % Turn on shadowing
    cfgModel.SampleRate=fs;
    cfgLayout = createUsersLayout(x(1),y(1),scenario);
        winChannel = comm.WINNER2Channel(cfgModel,cfgLayout);
        
        
        CSI = winner2.wim(cfgModel,cfgLayout);
        %
        CSI_Conv = cell2mat(CSI);
        
        
        Y=zeros(1,24);
        maxNumErrors = 5;   % The maximum number of packet errors at an SNR point
        maxNumPackets = 50; % The maximum number of packets at an SNR point
        
        packetErrorRate = zeros(1,3);
       
        
         iter = [1;2;3];
        tbl = table(iter);
        counter=1;
        startMcs =0;
        endMcs = 0;
        once = true;
         for mcs_ = startMcs:endMcs
            
             %changing MCS
            cfgHE.MCS=mcs_;
            for GI = [3.2,0.8]
                cfgHE.GuardInterval = GI;
            
        
             
            for isnr = 1:3 % iterations
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
            release(winChannel);
            end
            %% Plot Packet Error Rate vs SNR
        %     tbl = table(snr,Output1,Output2);
            MCS=packetErrorRate';
            tbl = addvars(tbl,MCS);
            % calculate average packet error rate above all snr
            % compare with threshold 10%
            %
            if mean(MCS,'all') < 0.1
        
                Y(counter)=1;
            %
            else 
                Y(counter)=0;
            end
        
            counter=counter+1;
            end
             
            
         end   
        disp(Y);
        data={abs(CSI_Conv),Y};
        f=data;
end