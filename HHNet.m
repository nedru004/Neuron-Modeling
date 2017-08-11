function  [spk_dark,spk_light]= HHNet(t, N, p)
if nargin<1
    t=200;
end;
if nargin<2
    N=1000;
end;
if nargin<3
    p=0.01;
end;

% Equations for activation and inactivation
alphaN = @(V) 0.01.*(V+55)./(1-exp(-(V+55)/10));
betaN = @(V) 0.125.*exp(-(V+65)/80);
alphaM = @(V)0.1.*(V+40)./(1-exp(-(V+40)/10));
betaM = @(V)4.*exp(-(V+65)/18.0);
alphaH = @(V)0.07.*exp(-(V+65)/20.0);
betaH = @(V) 1./(1+exp(-(V+35)/10));

dt = 0.01;
NTime = t/dt;

%define constants
gNa = 260;
gK = ones(N,1)*65;
%CHANGE CONDUCTANCE HERE
change = 0;
gK = gK+change;

gL = 0.6;
eNa = 50;     eK = -77;     eL = -55;
inhib = 0.2; % percentage of inhibitory neurons
V = -65*ones(N,1);
m=0.0532*ones(N, 1);
n=0.3183*ones(N, 1);
h=0.5947*ones(N, 1);
s=zeros(N, 1); %slow variable
f=zeros(N, 1); %fast variable

iSyn = -85;
eSyn = 0;   %reversal potnetial of the synapses 0=excitatory -85=inhibitory
taus =4;    %slow time constant of synapses
tauf = 1;   %fast time constant of synapses
gsyn= .1;   %global synaptic strength constant (individual weights are set below)

Input = 6; % current drive to all the cells
NoiseInputStrength = 6; % scales random independent current to individual cells


%% Generate random x and y positions for the neurons on a square 1x1
xpos = rand(N,1);
ypos = rand(N,1);

D=zeros(N);
for i=1:N
    for j=i+1:N
        D(i,j)=sqrt((xpos(i)-xpos(j))^2+(ypos(i)-ypos(j))^2); %calculate the Euclidian distance between the cells.
        D(j,i)=D(i,j);
    end;
end;

%% generate network connectivity
Nsyn = poissrnd(N*p,N, 1);  %choosing the number of synapses per cell from a poisson distirbution.

alpha = 0.5; %This is the space constant of probability of connecvity with distance.
meanSynStrength=1; %The mean synaptic strength
STDSynStrength=.1; %The standard deviation around the synaptic strength.
for i=1:N
    %now go select synapses.
    npost = 0;
    Syn{i}=zeros(1,Nsyn(i));                        % >>>>>>>>>>>>>>>>Understand how many synapse each cell has, how generated, average number of synapses
    SynWeights{i} = zeros(1,Nsyn(i));
    while(npost<Nsyn(i)) %While the number of selected post synaptic neurons is fewer than the number the cell is supposed to have, keep selecting synapses.
        p=0;
        r = 1;
        while(r>p)                     %keep randomly selecting cells and drawing a random number to determine if they are connnected until you find a connected pair.
            postcell = randperm(N,1);   % first, select a cell at random to decide if they are connected
            if(postcell~=i)             %if that cell is not the pre-synaptic cell, now we decide whether they are connected.
                r=rand(1);              %Draw a random number and decide if they are connected.
                p=exp(-D(i,postcell)/alpha); %Calculate the probability that they will be connected given their distance.
                if(r<=p)                %if r<p then we connect them
                    npost=npost+1;      %iterate the number of synapses this cell has.
                    Syn{i}(npost)= postcell;
                    SynWights{i}(npost)=randn(1)*STDSynStrength+meanSynStrength; %Select a synatpic strength at random
                    while(SynWights{i}(npost)<0) %Check that the synaptic strength is not negative, if it is, keep drawing strengths.
                        SynWights{i}(npost)=randn(1)*STDSynStrength+meanSynStrength;
                    end;
                end;
            end;                                            % >>>>>>>>>>>>>>>>>>autapses?
        end;
    end;
end;


%% run simulation
spk = zeros(1000,2);
spkcnt=0;
neuronType = ones(N, 1)*eSyn;              % create N-by-1 matrix to label cells as excitatory or inhibitory
neuronType(randi(N,floor(N*inhib))) = iSyn;                % where a=0 will represent an excitatory cell

for i=2:NTime  % Run time simulation
    VLast = V;
    % Hodgkin-Huxley Differential Equations
    dmdt = alphaM(V).*(1-m)-betaM(V).*m;                         %dm/dt
    dndt = alphaN(V).*(1-n)-betaN(V).*n;                         %dn/dt
    dhdt = alphaH(V).*(1-h)-betaH(V).*h;                         %dh/dt
    dsdt = -s/taus;
    dfdt = -f/tauf;
    
    if(i == 2000)
        pick = randi(N);
        s(pick)=s(pick)+1;
        f(pick)=f(pick)+1;
    end;
    dVdt = -gNa.*m.^3.*h.*(V-eNa)-gK.*n.^4.*(V-eK)-gL.*(V-eL)+gsyn.*(s-f).*(neuronType-V)+Input +NoiseInputStrength*randn(N,1);
    
    
    % Euler integrator
    V = V+dt.*dVdt;
    m = m+dt.*dmdt;
    n = n+dt.*dndt;
    h = h+dt.*dhdt;
    s = s+dt.*dsdt;
    f = f+dt.*dfdt;
    
    spkind = find((VLast<0).*(V>0)); %% list indeces of cells whose voltage was below 0 volts last time step and above this time step
    for cind=1:length(spkind)
        spkcnt = spkcnt+1;
        if(spkcnt>length(spk))
            spk = [spk; zeros(1000,2)];
        end;
        
        spk(spkcnt,1)=i*dt; % Multiply by dt to get spike time in milliseconds
        spk(spkcnt,2)=spkind(cind); % Record which cell spiked
        
        if ~isempty(Syn{cind}) % This strengthens synaptic connections that are firing
            s(Syn{cind}) = s(Syn{cind}) +SynWights{cind}';  % update the slow variables of the postsynatpic neurons
            f(Syn{cind}) = f(Syn{cind}) +SynWights{cind}';  % update the fast variables of the postsynatpic neurons
        end;
    end;
end;

spk_dark = spk;
% Simulate Blue Light Stimulation
change = lognrnd(1,.75,N,1); % This simulates a Lumitoxin
gK = gK+change;

disp('Finished Dark Simulation')

%% run second simulation
spk = zeros(1000,2);
spkcnt=0;
for i=2:NTime  % Run time simulation
    VLast = V;
    %Hodgkin-Huxley Differential Equations
    dmdt = alphaM(V).*(1-m)-betaM(V).*m;                         %dm/dt
    dndt = alphaN(V).*(1-n)-betaN(V).*n;                         %dn/dt
    dhdt = alphaH(V).*(1-h)-betaH(V).*h;                         %dh/dt
    dsdt = -s/taus;
    dfdt = -f/tauf;
    
    if(i == 2000)
        pick = randi(N);
        s(pick)=s(pick)+1;
        f(pick)=f(pick)+1;
    end;
    dVdt = -gNa.*m.^3.*h.*(V-eNa)-gK.*n.^4.*(V-eK)-gL.*(V-eL)+gsyn.*(s-f).*(neuronType-V)+Input +NoiseInputStrength*randn(N,1);
    
    
    %Euler integrator
    V = V+dt.*dVdt;
    m = m+dt.*dmdt;
    n = n+dt.*dndt;
    h = h+dt.*dhdt;
    s = s+dt.*dsdt;
    f = f+dt.*dfdt;
    
    spkind = find((VLast<0).*(V>0)); % list indeces of cells whose voltage was below 0 volts last time step and above this time step ::: AKA neuron spike
    for cind=1:length(spkind)
        spkcnt = spkcnt+1;
        if(spkcnt>length(spk))
            spk = [spk; zeros(1000,2)];
        end;
        spk(spkcnt,1)=i*dt;
        spk(spkcnt,2)=spkind(cind);
        
        if(length(Syn{cind}))
            s(Syn{cind}) = s(Syn{cind}) +SynWights{cind}';  %% update the slow variables of the postsynatpic neurons
            f(Syn{cind}) = f(Syn{cind}) +SynWights{cind}';  %% update the fast variables of the postsynatpic neurons
        end;
    end;
end;

disp('Finished Light Simulation')

spk_light = spk;

%% plot results
% spk = spk(1:spkcnt,:);
% plot(dt*spk(:,1),spk(:,2),'.');
% xlabel('Time (msec)');
% ylabel('Cell number');
% title('Spike rastergram');
% spkhist = hist(spk(:,1), t); % calculates the spike densities on a millisecond time scale.
% spkhist = spkhist(floor(t*3/4):end);
% phist=spkhist/sum(spkhist);
% ind = phist>0;
% ent = sum(phist(ind).*log(phist(ind)));
% EEnt=log(1/t);
% Nent = (EEnt-ent)/EEnt; % Gives a normalized entropy where 1 means synchronized and zero means uniformly distributed.


save([date,'_HHNet.mat'],'spk_dark','spk_light','change','neuronType')
