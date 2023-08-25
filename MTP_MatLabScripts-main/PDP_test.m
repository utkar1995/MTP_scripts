figure(Position=[scrsz(3)*.3-figSize/2,scrsz(4)*.25-figSize/2,figSize,figSize]);
        hold on;
        for linkIdx = 1:1
            delay = chanInfo.ChannelFilterDelay(linkIdx);
            stem(((0:(frameLen-1))-delay)/chanInfo.SampleRate(linkIdx), ...
                abs(rx{linkIdx}(1:1600,1)));
        end
        maxX = max((cell2mat(cellfun(@(x) find(abs(x) < 1e-8,1,"first"), ...
            rx.',UniformOutput=false)) - chanInfo.ChannelFilterDelay)./ ...
            chanInfo.SampleRate);
        minX = -max(chanInfo.ChannelFilterDelay./chanInfo.SampleRate);
        xlim([minX, maxX]);
        xlabel("Time (s)"); 
        ylabel("Magnitude");
        legend("Link 1");
        title("Impulse Response at First Receive Antenna");