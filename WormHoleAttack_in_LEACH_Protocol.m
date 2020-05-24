close all;
clear all
clear all;
clc
%% Initial parameters
Area_x=200; % Required Area Dimensions
Area_y=200;

no_of_sensors=50;% Sensors amount
n=no_of_sensors;
SinkNode_x=100; %coordinates of the SinkNode
SinkNode_y=190;

x=0; y=0; %temp variable for better result
dead_nodes=0; % Initally the number of death node will be zero
             %Energy Requirement(units in Joules/bit)
Eo=2;
Eelec=50*10^(-9);
ETx=50*10^(-9);
ERx=50*10^(-9);
 
%Transmit Amplifier Types
Eamp=100*10^(-12);
%Data Aggregation Energy
EDA=5*10^(-9);
%Size of the data package
k=4000; 
%
%let add 5 millicious nodes in the network
no_of_malicious_node=2;
%let 5% of the total amount of nodes used in the network is proposed to
%give good resuls
p=0.05;
%Number of Cluster
No=p*n;
%Round of Operation
rnd=0;
%Details of current num.of.operating nodes
operating_nodes=n;
transmissions=0;
temp_val=0;
flag1stdead=0;
tcount=1;

%% Creating WSensorNode Network of normal nodes

for i=1:no_of_sensors
    SensorNode(i).id=i;
    SensorNode(i).x=rand(1,1)*Area_x;
    SensorNode(i).y=rand(1,1)*Area_y;
    SensorNode(i).E=Eo; %initial Energy of each sensor Node
    SensorNode(i).malicious='n';
    SensorNode(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SensorNode(i).cluster=0;	% the cluster which a node belongs to
    SensorNode(i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SensorNode(i).rop=0;	% number of rounds node was operational
    SensorNode(i).rleft=0;  % rounds left for node to become available for Cluster Head election
    SensorNode(i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
    SensorNode(i).dts=0;    % nodes distance from the sink
    SensorNode(i).tel=0;	% states how many times the node was elected as a Cluster Head
    SensorNode(i).rn=0;     % round node got elected as cluster head
    SensorNode(i).chid=0;
    hold on;
    figure(1);
    plot(x,y,Area_x,Area_y,SensorNode(i).x,SensorNode(i).y,'ob',SinkNode_x,SinkNode_y,'*r');
    title ({'WormHole Attack on LEACH'; 'Wireless Sensor network';})
    xlabel 'x-axis(in meter)';
    ylabel 'y-axis(in meter)';
    grid on;
end

%% Add malicious nodes into the network
% let make node 10,20 as a millicious node
for i=1:no_of_malicious_node
    SensorNode(i*10).malicious='m';
    plot(x,y,Area_x,Area_y,SensorNode(i*10).x,SensorNode(i*10).y,'or');
    hold on;
end
mx=[SensorNode(10).x SensorNode(20).x];
my=[SensorNode(10).y SensorNode(20).y];
plot(mx,my);
hold on;
text(mx(1),my(1)-0.05,'Worm Hole tunnel','HorizontalAlignment','center');


%% operating phase

while operating_nodes>0        
    % Displays Current Round %     
    rnd     
	% Threshold Value %
	t=(p/(1-p*(mod(rnd,1/p))));
    % Re-election Value %
    tleft=mod(rnd,1/p);
	% Reseting Previous Amount Of Cluster Heads In the Network %
	CLheads=0;
    % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
    energy=0;        
% Cluster Heads Election %
        for i=1:n
            SensorNode(i).cluster=0;    % reseting cluster in which the node belongs to
            SensorNode(i).role=0;       % reseting node role
            SensorNode(i).chid=0;       % reseting cluster head id
            if SensorNode(i).rleft>0
               SensorNode(i).rleft=SensorNode(i).rleft-1;
            end
            if (SensorNode(i).E>0) && (SensorNode(i).rleft==0)
                generate=rand;	
                    if generate< t  %Threshold Value(t)
                    SensorNode(i).role=1;	% assigns the node role of acluster head
                    SensorNode(10).clusteridinfo(tcount).data=SensorNode(i).id;
                    SensorNode(20).clusteridinfo(tcount).data=SensorNode(i).id;
                    tcount=tcount+1;
                    SensorNode(i).rn=rnd;	% Assigns the round that the cluster head was elected to the data table 
                     SensorNode(i).tel=SensorNode(i).tel + 1;   
                    SensorNode(i).rleft=1/p-tleft;    % rounds for which the node will be unable to become a CH
                    SensorNode(i).dts=sqrt((SinkNode_x-SensorNode(i).x)^2 + (SinkNode_y-SensorNode(i).y)^2); % calculates the distance between the sink and the cluster hea
                    CLheads=CLheads+1;	% sum of cluster heads that have been elected 
                    SensorNode(i).cluster=CLheads; % cluster of which the node got elected to be cluster head
                    CL(CLheads).x=SensorNode(i).x; % X-axis coordinates of elected cluster head
                    CL(CLheads).y=SensorNode(i).y; % Y-axis coordinates of elected cluster head
                    CL(CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
                    end
            end
        end
        %Totalnumber_of_cluster=CLheads;
	% Fixing the size of "CL" array %
	CL=CL(1:CLheads);
% Grouping the Nodes into Clusters & caclulating the distance between node and cluster head %
     CLheads
       for i=1:n
        if  (SensorNode(i).role==0) && (SensorNode(i).E>0) && (CLheads>0) % if node is normal
            for m=1:CLheads
            d(m)=sqrt((CL(m).x-SensorNode(i).x)^2 + (CL(m).y-SensorNode(i).y)^2);
            % we calculate the distance 'd' between the sensor node that is
            % transmitting and the cluster head that is receiving with the following equation+ 
            % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
            % the cluster head and x1 and y1 the coordinates of the transmitting node
            end
        d=d(1:CLheads); % fixing the size of "d" array
        [M,I]=min(d(:)); % finds the minimum distance of node to CH
        [Row, Col] = ind2sub(size(d),I);
        % displays the Cluster Number in which this node belongs too
        SensorNode(i).cluster=Col; % assigns node to the cluster
        SensorNode(i).dtch= d(Col); % assigns the distance of node to CH
        SensorNode(i).chid=CL(Col).id;
        end
       end
                           %%%%%% Steady-State Phase %%%%%%
% Energy Dissipation for normal nodes %    
    for i=1:n
       if (SensorNode(i).cond==1) && (SensorNode(i).role==0) && (CLheads>0)
       	if SensorNode(i).E>0
            ETx= Eelec*k + Eamp * k * SensorNode(i).dtch^2;
            SensorNode(i).E=SensorNode(i).E - ETx;
            energy=energy+ETx;
        % Dissipation for cluster head during reception
        if SensorNode(SensorNode(i).chid).E>0 && SensorNode(SensorNode(i).chid).cond==1 && SensorNode(SensorNode(i).chid).role==1
            ERx=(Eelec+EDA)*k;
            energy=energy+ERx;
            SensorNode(SensorNode(i).chid).E=SensorNode(SensorNode(i).chid).E - ERx;
             if SensorNode(SensorNode(i).chid).E<=0  % if cluster heads energy depletes with reception
                SensorNode(SensorNode(i).chid).cond=0;
                SensorNode(SensorNode(i).chid).rop=rnd;
                dead_nodes=dead_nodes +1;
                operating_nodes= operating_nodes - 1
             end
        end
        end
        
        
        if SensorNode(i).E<=0       % if nodes energy depletes with transmission
        dead_nodes=dead_nodes +1;
        operating_nodes= operating_nodes - 1
        SensorNode(i).cond=0;
        SensorNode(i).chid=0;
        SensorNode(i).rop=rnd;
        end
        
      end
    end            
    
    
    
% Energy Dissipation for cluster head nodes %
   
   for i=1:n
     if (SensorNode(i).cond==1)  && (SensorNode(i).role==1)
         if SensorNode(i).E>0
            ETx= (Eelec+EDA)*k + Eamp * k * SensorNode(i).dts^2;
            SensorNode(i).E=SensorNode(i).E - ETx;
            energy=energy+ETx;
         end
         if  SensorNode(i).E<=0     % if cluster heads energy depletes with transmission
         dead_nodes=dead_nodes +1;
         operating_nodes= operating_nodes - 1
         SensorNode(i).cond=0;
         SensorNode(i).rop=rnd;
         end
     end
   end

   

  
    if operating_nodes<n && temp_val==0
        temp_val=1;
        flag1stdead=rnd
    end
    % Display Number of Cluster Heads of this round %
    CLheads;
   
    
    transmissions=transmissions+1;
    if CLheads==0
    transmissions=transmissions-1;
    end
    
 
    % Next Round %
    rnd= rnd +1;
    
    tr(transmissions)=operating_nodes;
    op(rnd)=operating_nodes;
    

    if energy>0
    nrg(transmissions)=energy;
    end
    

end


sum=0;
for i=1:flag1stdead
    sum=nrg(i) + sum;
end

temp1=sum/flag1stdead;
temp2=temp1/n;

for i=1:flag1stdead
avg_node(i)=temp2;
end
%% plotting simulation Results    
    % Plotting Simulation Results "Operating Nodes per Round" %
    figure(2)
    plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
    title ({'WormHole Attack on LEACH'; 'Operating Nodes per Round';})
    xlabel 'Rounds';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(3)
    plot(1:transmissions,tr(1:transmissions),'-r','Linewidth',2);
    title ({'WormHole Attack on LEACH'; 'Operational Nodes per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(4)
    plot(1:flag1stdead,nrg(1:flag1stdead),'-r','Linewidth',2);
    title ({'WormHole Attack on LEACH'; 'Energy consumed per Transmission';})
    xlabel 'Transmission';
    ylabel 'Energy ( J )';
    hold on;
    

    % Plotting Simulation Results  %
    figure(5)
    plot(1:flag1stdead,avg_node(1:flag1stdead),'-r','Linewidth',2);
    title ({'WormHole Attack on LEACH'; 'Average Energy consumed by a Node per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Energy ( J )';
    hold on;
  
  

