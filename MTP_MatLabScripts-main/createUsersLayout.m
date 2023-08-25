function f = createUsersLayout(x,y,scenario,antennas)
    rng(200);
    users=1;

    
    BSAA = winner2.AntennaArray("UCA",antennas,0.02);   % UCA-8 antenna array for BS
    arrayOfNames=[BSAA];
    
    for i= 1:users
        arrayOfNames=[arrayOfNames,winner2.AntennaArray("UCA",antennas,0.05)];
    end
    MSIdx = linspace(2,users+1,users); 
    BSIdx = {1};
    NL = antennas;
    maxRange = 100;
    rndSeed = 101;
    
    cfgLayout = winner2.layoutparset(MSIdx,BSIdx,NL, ...
       arrayOfNames,maxRange,rndSeed);
%     x=randi([1 170],1,users);
%     y=randi([1 100],1,users);
    cfgLayout.Stations(1).Pos(1:2) = [100, 100];
    for i = 1:users
        cfgLayout.Stations(2).Pos(1:2) = [x, y];
        v = rand(3,1) - 0.5;
        cfgLayout.Stations(i+1).Velocity = v/norm(v,'fro');
    end
    
    
    row1=ones(1,users);
    row2 = linspace(2,users+1,users); 
    cfgLayout.Pairing = [row1;row2];
    cfgLayout.PropagConditionVector = ones(1,users);
    cfgLayout.ScenarioVector = [scenario];
    f=cfgLayout;

end

