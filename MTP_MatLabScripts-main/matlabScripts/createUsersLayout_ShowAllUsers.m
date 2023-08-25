function f = createUsersLayout_ShowAllUsers(users,scenario)
    rng(200);
    

    
    BSAA = winner2.AntennaArray("UCA",1,0.02);   % UCA-8 antenna array for BS
    arrayOfNames=[BSAA];
    
    for i= 1:users
        arrayOfNames=[arrayOfNames,winner2.AntennaArray("ULA",1,0.005)];
    end
    MSIdx = linspace(2,users+1,users); 
    BSIdx = {1};
    NL = 2;
    maxRange = 100;
    rndSeed = 101;
    
    cfgLayout = winner2.layoutparset(MSIdx,BSIdx,NL, ...
       arrayOfNames,maxRange,rndSeed);
    x=randi([20 300],1,users);
    y=randi([20 300],1,users);
    cfgLayout.Stations(1).Pos(1:2) = [150, 150];
    for i = 1:users
        cfgLayout.Stations(i+1).Pos(1:2) = [x(i), y(i)];
    end
    
    
    row1=ones(1,users);
    row2 = linspace(2,users+1,users); 
    cfgLayout.Pairing = [row1;row2];
    cfgLayout.PropagConditionVector = zeros(1,users);
    cfgLayout.ScenarioVector = [scenario];
    f=cfgLayout;
    numBSSect = sum(cfgLayout.NofSect);
    %visualizing MS BS positions
    %
    BSPos = cell2mat({cfgLayout.Stations(1:numBSSect).Pos});
    MSPos = cell2mat({cfgLayout.Stations(numBSSect+1:end).Pos});
    
    scrsz = get(groot,"ScreenSize");
    figSize = min(scrsz([3,4]))/2.3;
    figure( ...
        Position=[scrsz(3)*.5-figSize/2,scrsz(4)*.7-figSize/2,figSize,figSize]);
    hold on; 
    grid on;
    hBS = plot(BSPos(1,:),BSPos(2,:),"or");   % Plot BS
    hMS = plot(MSPos(1,:),MSPos(2,:),"xb");   % Plot MS
    for linkIdx = 1:users                 % Plot links
        pairStn = cfgLayout.Pairing(:,linkIdx);
        pairPos = cell2mat({cfgLayout.Stations(pairStn).Pos});
        plot(pairPos(1,:),pairPos(2,:),"-b");
    end
    xlim([0 300]); ylim([0 300]);
    xlabel("X Position (meters)");
    ylabel("Y Position (meters)")
    legend([hBS, hMS],"BS","MS",location="northwest");
    % end visualization


end

