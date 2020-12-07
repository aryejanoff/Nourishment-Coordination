function [R2coord_star,R3coord_star,R2uncoord_star_cons,R3uncoord_star_cons,R2uncoord_star_risk,R3uncoord_star_risk,M_coord,M_uncoord_cons,M_uncoord_risk,NB2_coord,NB3_coord,NB2_uncoord_cons,NB3_uncoord_cons,NB2_uncoord_risk,NB3_uncoord_risk,behavior_coord,behavior_uncoord_cons,behavior_uncoord_risk,Neff_coord,Neff_uncoord_cons,Neff_uncoord_risk]=maincode(alpha2,alpha3)

%Nourishment Rotation Length Outer Loops
R_min=0.2; %minimum rotation length (yrs)
R_max=25; %maximum rotation length (yrs)
dR=0.2; %vector step between rotation lengths
R2_vec=R_min:dR:R_max; R3_vec=R_min:dR:R_max; %rotation length vectors (only nourishment options)
R2_vector=[NaN R2_vec(1:end)]; R3_vector=[NaN R3_vec(1:end)]; %rotation length vectors (nourish/no-nourish options)
n2=length(R2_vector); n3=length(R3_vector); %length of nourishment vectors
NB2_storage=NaN(n2,n3); NB3_storage=NaN(n2,n3); TNB_storage=NaN(n2,n3); %community-specific and system net benefit matrix storage
behavior_storage=NaN(n2,n3); %emergent behaviors matrix storage (i.e., seaward growth, hold the line, slow retreat, full retreat)
N_efficiency_storage=NaN(n2,n3); %physical nourishment efficiency matrix storage

parfor iR2=1:numel(R2_vector) %loop through possible management options in one community for optimization
    NB2_vector_R=NaN(1,n3); NB3_vector_R=NaN(1,n3); TNB_vector_R=NaN(1,n3); %community-specific and system net benefit vector storage (suggested for use with parfor) 
    behavior_vector=NaN(1,n3); %emergent behaviors vector storage (suggested for use with parfor) 
    N_efficiency_vector=NaN(1,n3); %physical nourishment efficiency vector storage (suggested for use with parfor) 
    
    for iR3=1:numel(R3_vector) %loop through possible management options in other community for optimization 
        %% Input Physical Pmarameters %%
        lot_size=30; %length of each property assuming a square property (m)
        w_init=30; %initial beach width (m)
        beta=0.4; %beach width hedonic parameter 
%         alpha2=5e6/(w_init^beta); %baseline property value in community 1, i.e., value less the environmental beach amenity ($)
%         alpha3=5e6/(w_init^beta); %baseline property value in community 2
        s=1500; %alongshore community length (m)
        rows_along=s/lot_size; %number of beachfront properties in a community alongshore
        rows_cross=1; %number of cross-shore property rows
        properties_total=rows_cross*rows_along; %total number of properties in a community
        comm_width=rows_cross*lot_size; %initial community width (m)
        psi=0.2; %hedonic parameter for cross-shore property value change as distance from ocean increases (unused in paper but built in as an extension for future analysis)
        D=16; %depth of closure (m)
        gamma=5; %background erosion rate (m/yr)
        d=500000; %alongshore flux coefficient (m^2/yr)
        K=2000; %cross-shore flux coefficient (m^2/yr)
        phi=10; %sand unit cost ($/m^3) 
        c=2e6; %fixed nourishment cost ($)
        rho=0.06; %discount rate (yr^-1)
        R2=R2_vector(iR2); %rotation length for community 1
        R3=R3_vector(iR3); %rotation length for community 2
        xN=50; %nourishment magnitude (m)
        nu=0; %beach width decline beyond max threshold (unused in this paper but built in as an extension for future analysis)
        theta_eq=0.02; %equilibrium shoreface slope (m/m)
        k2=0; %nourishment counter in community 1
        k3=0; %nourishment counter in community 2
        
        %% Computational Parameters %%
        tmax=50; dt=0.05; t=0:dt:tmax; n=length(t); %time horizon (yrs); time step (yrs); time vector; length of time vector
        Smax=3; ds=1; S=1:ds:Smax; m=length(S); %number of alongshore cells in model domain; space step; space vector; length of space vector
        A2=alpha2*rho; A3=alpha3*rho; %baseline rental value for each community using simplified amortization formula ($/yr)
        theta=zeros(n,m); qL=zeros(n,m); qC=zeros(n,m); fS=zeros(n,m); fT=zeros(n,m); %matrix storage for: shoreface slope (m/m); alongshore flux (m^2/yr); cross-shore flux (m^2/yr); shoreline change ODE (m/yr); shoreface toe change ODE (m/yr)
        xS=zeros(n,m); xT=zeros(n,m); xH=zeros(n,m); w=zeros(n,m); %storage: shoreline location (m); shoreface toe location (m); property setback (m); beach width (m)
        B=zeros(n,m); B_disc=zeros(n,m); TB=zeros(n,m); %storage: continuous benefits ($/yr); discounted benefits ($); total benefits ($)
        C=zeros(n,m); TC=zeros(n,m); %storage: discounted discrete costs ($); total costs ($)
        nb=NaN(n,m); N_h=zeros(n,m); %storage: discounted and integrated net benefits ($); number of cross-shore property rows through time       
        volume=zeros(n,m); nE=zeros(n,m); %time-specific nourishment volumes (m^3); tracker for number of nourishment episodes
        q_loss=zeros(1,n); q_along=zeros(1,n); q_cross=zeros(1,n); q_gamma=zeros(1,n); %change rates for total sediment loss; alongshore sediment loss; cross-shore sediment loss; and loss due to background erosion (m/yr)
        V_nourish=zeros(1,n); V_loss=zeros(1,n); efficiency=zeros(1,n); %total nourishment volume across both communities (m^3); total sand lost from both communites due to cross-shore/alongshore dynamics and background erosion (m^3); system-wide nourishment efficiency 
        behavior=NaN(1); NB2=NaN(1); NB3=NaN(1); TNB=NaN(1); %create variables for storing in external scripts: emergent behavior; community 1 net benefit ($); community 2 net benefit ($); total net benefit for both communities ($)

        %% Initial Conditions %%
        xS(1,:)=comm_width+w_init; %initial shoreline location assuming equilibrium shoreface slope and initial beach width (m)
        xT(1,:)=xS(1,:)+(D/theta_eq); %initial shoreface toe location assuming equilibrium shoreface slope and initial beach width (m)
        xH(1,:)=comm_width; %initial property setback location (m)
        N_h(1,:)=rows_cross; %initial number of cross-shore property rows
        w(1,:)=xS(1,:)-xH(1,:); %initial beach width (m)

        %% Main Code %%
        for i=1:n-1 %time for loop
            for j=2:m-1 %space for loop                          
                %% Nourishment Initiation + Volume (see input and computational parameters sections for variable/parameter descriptions)
                ff2=round(t(i)-k2*R2,4);
                if ff2==0 %populates nourishment volume matrix for community 1 based on its rotation length
                    k2=k2+1; volume(i,2)=0.5*s*D*xN; nE(i,2)=1;
                end
                ff3=round(t(i)-k3*R3,4);
                if ff3==0 %populates nourishment volume matrix for community 2 based on its rotation length
                    k3=k3+1; volume(i,3)=0.5*s*D*xN; nE(i,3)=1;
                end
                V_nourish(i+1)=V_nourish(i)+(2*(volume(i,2)+volume(i,3))/(s*D));
                
                %% Fluxes (Along/Cross-shore) and Shoreface Dynamics (see input and computational parameters sections for variable/parameter descriptions)
                qL(i,j)=d*((xS(i,j-1)-xS(i,j))/s); qL(i,1)=d*((xS(i,m)-xS(i,1))/s); qL(i,m)=d*((xS(i,m-1)-xS(i,m))/s); 
                theta(i,j)=D/(xT(i,j)-xS(i,j)); theta(i,1)=D/(xT(i,1)-xS(i,1)); theta(i,m)=D/(xT(i,m)-xS(i,m)); 
                qC(i,j)=K*(theta(i,j)-theta_eq); qC(i,1)=K*(theta(i,1)-theta_eq); qC(i,m)=K*(theta(i,m)-theta_eq);
                q_along(i)=((2*(qL(i,1)-qL(i,2)))/(s(2)+s(3))); q_cross(i)=((4*(qC(i,2)+qC(i,3)))/D); q_gamma(i)=gamma;
                q_loss(i)=q_along(i)+q_cross(i)+gamma; 
                V_loss(i+1)=V_loss(i)+dt*q_loss(i);
                efficiency(i)=V_nourish(i)./(V_nourish(i)+V_loss(i));
                
                %% ODE's (see input and computational parameters sections for variable/parameter descriptions)
                fT(i,j)=(4*K*(theta(i,j)-theta_eq)/D)-gamma; fT(i,1)=(4*K*(theta(i,1)-theta_eq)/D)-gamma; fT(i,m)=(4*K*(theta(i,m)-theta_eq)/D)-gamma;
                fS(i,j)=((2*d/(s^2))*(xS(i,j-1)-2*xS(i,j)+xS(i,j+1)))-(4*K*(theta(i,j)-theta_eq)/D)-gamma; fS(i,1)=((2*d/(s^2))*(xS(i,m)-2*xS(i,1)+xS(i,2)))-(4*K*(theta(i,1)-theta_eq)/D)-gamma; fS(i,m)=((2*d/(s^2))*(xS(i,m-1)-2*xS(i,m)+xS(i,1)))-(4*K*(theta(i,m)-theta_eq)/D)-gamma; 

                %% Numerical Approximations using Forward Euler Method (see input and computational parameters section for variable/parameter descriptions)

                xT(i+1,j)=xT(i,j)+dt*fT(i,j); xT(i+1,1)=xT(i,1)+dt*fT(i,1); xT(i+1,m)=xT(i,m)+dt*fT(i,m);
                xS(i+1,1)=xS(i,1)+dt*fS(i,1); 
                if xS(i,j)<=lot_size %stops nourishment if property is lost in cell j
                    volume(i,j)=0; volume(i+1,j)=0; N_h(i,j)=0; xH(i,j)=0; N_h(i+1,j)=0; xH(i+1,j)=0;  
                else
                    xH(i+1,j)=xH(i,j);
                end
                if xS(i,m)<=lot_size %stops nourishment if property is lost in cell m
                    volume(i,m)=0; volume(i+1,m)=0; N_h(i,m)=0; xH(i,m)=0; N_h(i+1,m)=0; xH(i+1,m)=0;  
                else
                    xH(i+1,m)=xH(i,m);
                end
                if volume(i,2)~=0 %adds nourishment sand to shoreline based on volume matrix for community 1
                    xS(i+1,2)=xS(i,2)+xN;
                elseif volume(i,2)==0
                    xS(i+1,2)=xS(i,2)+dt*fS(i,2);
                end
                if volume(i,3)~=0 %adds nourishment sand to shoreline based on volume matrix for community 2
                    xS(i+1,3)=xS(i,3)+xN; 
                elseif volume(i,3)==0
                    xS(i+1,3)=xS(i,3)+dt*fS(i,3);
                end
                w(i,j)=xS(i,j)-xH(i,j); w(i,1)=xS(i,1)-xH(i,1); w(i,m)=xS(i,m)-xH(i,m);

                %% Housing Lines 
                N_h(i,j)=xH(i,j)/lot_size; N_h(i,1)=xH(i,1)/lot_size; N_h(i,m)=xH(i,m)/lot_size;
                
                %% Benefit 
                B(i,2)=rows_along*A2*((w(i,2)).^beta)*((N_h(i,2)).^psi)-(nu*(w(i,2).^2)); B(i,3)=rows_along*A3*((w(i,3)).^beta)*((N_h(i,3)).^psi)-(nu*(w(i,3).^2));
                B_disc(i,j)=B(i,j)*exp(-rho*t(i)); B_disc(i,m)=B(i,m)*exp(-rho*t(i));
                TB(i,j)=(dt/3)*(B_disc(1,j)+4*(sum(B_disc(2:2:end,j)))+2*sum(B_disc(2:1:end,j))+B_disc(end,j)); TB(i,m)=(dt/3)*(B_disc(1,m)+4*(sum(B_disc(2:2:end,m)))+2*sum(B_disc(2:1:end,m))+B_disc(end,m));

                %% Cost 
                C(i,j)=(c+volume(i,j)*phi)*exp(-rho*t(i)); C(i,m)=(c+volume(i,m)*phi)*exp(-rho*t(i));
                if volume(i,j)==0 %ensures no cost is incurred when there is no nourishment in cell j
                    C(i,j)=0;
                end
                if volume(i,m)==0 %ensures no cost is incurred when there is no nourishment in cell m
                    C(i,m)=0;
                end
                if N_h(i,j)==0 %ensures no cost is incurred when there are no properties in cell j
                    C(i,j)=0;
                end
                if N_h(i,m)==0 %ensures no cost is incurred when there are no properties in cell m
                    C(i,m)=0;
                end
                TC(i,j)=sum(C(:,j)); TC(i,m)=sum(C(:,m));

                %% marginal net benefit 
                nb(i,j)=TB(i,j)-TC(i,j); nb(i,m)=TB(i,m)-TC(i,m);

                %% Net Benefit
                NB2=nb(i,2); NB3=nb(i,3); TNB=NB2+NB3;
            end
        end
        
        %% Emergent Behavior Identification
        if (max(xS((tmax/dt-50):end,2))>comm_width+w_init+xN && max(xS((tmax/dt-50):end,3))>comm_width+w_init+xN)
            behavior=0; %seaward growth
        elseif max(xS((tmax/dt-50):end,2))<=comm_width+w_init+xN && max(xS((tmax/dt-50):end,3))<=comm_width+w_init+xN && xH(end,2)-xH(1,2)==0 && xH(end,3)-xH(1,3)==0
            behavior=3; %hold the line
        elseif  (xH(end,2)-xH(1,2)<0 && xH(end,3)-xH(1,3)<0) && (isnan(R2)==0 || isnan(R3)==0)
            behavior=6; %slow retreat
        elseif xH(end,2)==0 && xH(end,3)==0 && isnan(R2)==1 && isnan(R3)==1 
            behavior=9; %retreat
        elseif (max(xS((tmax/dt-50):end,2))>comm_width+w_init+xN || max(xS((tmax/dt-50):end,3))>comm_width+w_init+xN) && ((max(xS((tmax/dt-50):end,2))<=comm_width+w_init+xN && xH(end,2)-xH(1,2)==0) || (max(xS((tmax/dt-50):end,3))<=comm_width+w_init+xN && xH(end,3)-xH(1,3)==0))
            behavior=1.5; %mixing: seaward growth/hold the line
        elseif ((max(xS((tmax/dt-50):end,2))<=comm_width+w_init+xN && xH(end,2)-xH(1,2)==0) || (max(xS((tmax/dt-50):end,3))<=comm_width+w_init+xN && xH(end,3)-xH(1,3)==0)) && (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0)
            behavior=4.5; %mixing: hold the line/slow retreat
        elseif (max(xS((tmax/dt-50):end,2))>comm_width+w_init+xN || max(xS((tmax/dt-50):end,3))>comm_width+w_init+xN) && (xH(end,2)-xH(1,2)<0 || xH(end,3)-xH(1,3)<0)
            behavior=10.5; %mixing: seaward growth/slow retreat
        end

        %% Nourishment Efficiency
        if efficiency(end-1)==0
            Eff=NaN;
        elseif efficiency(end-1)~=0
            Eff=efficiency(end-1);
        end


        %% # nourishment episodes
        TnE2=sum(nE(:,2));
        TnE3=sum(nE(:,3));

        %% Output Storage Vectors
        N_efficiency_vector(iR3)=Eff;
        NB2_vector_R(iR3)=NB2;
        NB3_vector_R(iR3)=NB3;
        TNB_vector_R(iR3)=TNB;
        behavior_vector(iR3)=behavior;
    end
    N_efficiency_storage(iR2,:)=N_efficiency_vector;
    NB2_storage(iR2,:)=NB2_vector_R;
    NB3_storage(iR2,:)=NB3_vector_R;
    TNB_storage(iR2,:)=TNB_vector_R;
    behavior_storage(iR2,:)=behavior_vector;
end

%% find optimal rotation lengths and corresponding behaviors, efficiencies, and net benefits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coordination scheme
M_coord=max(TNB_storage(:)); %max net benefit in matrix
[row_coord,col_coord]=find(TNB_storage==M_coord); %finds location in matrix where maximum net benefit occurs
row_c=row_coord(1); col_c=col_coord(1); %ensures only one combination is output to external script in case optimization doesn't converge due to resolution issue (only occurs under specific parameter combinations and resolution issue is confirmed as the cause)
R2coord_star=R2_vector(row_c); %optimal rotation length in community 1
R3coord_star=R3_vector(col_c); %optimal rotation length in community 2
NB2_coord=NB2_storage(row_c,col_c); %net benefit realized by community 1
NB3_coord=NB3_storage(row_c,col_c); %net benefit realized by community 2
behavior_coord=behavior_storage(row_c,col_c); %emergent behavior for combination of optimal rotation lengths under coordination
Neff_coord=N_efficiency_storage(row_c,col_c); %nourishment efficiency for combination of optimal rotation lengths under coordination

%non-coordination schemes
M2=max(NB2_storage(:)); %maximum community 1 net benefit assuming neighbor nourishes as much as possible (risky assumption)
[row2,col2]=find(NB2_storage==M2); %location in net benefit matrix where above maximum occurs
row2_unrisk=row2(1); col2_unrisk=col2(1); %ensures only one selection occurs in case there is a resolution issue
R2uncoord_star_risk=R2_vector(row2_unrisk); %optimal rotation length in community 1 under risky non-coordination
R2_uncoord_assump=R3_vector(col2_unrisk); %community 1 assumption of neighbor's choice (i.e., that they will nourish as much as possible)
M2_cons=max(NB2_storage(:,1)); %maximum community 1 net benefit assuming neighbor does nothing (conservative assumption)
[row2_cons,col2_cons]=find(NB2_storage==M2_cons); %location in net benefit matrix where above maximum occurs
row2_uncons=row2_cons(1); col2_uncons=col2_cons(1); %ensures only one selection occurs in case there is a resolution issue
R2uncoord_star_cons=R2_vector(row2_uncons); %optimal rotation length in community 1 under conservative non-coordination
R2_uncoord_assump_cons=R3_vector(col2_uncons); %community 1 assumption of neighbor's choice (i.e., that they will not nourish)

M3=max(NB3_storage(:)); %maximum community 2 net benefit assuming neighbor nourishes as much as possible (risky assumption)
[row3,col3]=find(NB3_storage==M3); %location in net benefit matrix where above maximum occurs
row3_unrisk=row3(1); col3_unrisk=col3(1); %ensures only one selection occurs in case there is a resolution issue
R3uncoord_star_risk=R3_vector(col3_unrisk); %optimal rotation length in community 2 under risky non-coordination
R3_uncoord_assump=R2_vector(row3_unrisk); %community 2 assumption of neighbor's choice (i.e., that they will nourish as much as possible)
M3_cons=max(NB3_storage(1,:)); %maximum community 2 net benefit assuming neighbor does nothing (conservative assumption)
[row3_cons,col3_cons]=find(NB3_storage==M3_cons); %location in net benefit matrix where above maximum occurs
row3_uncons=row3_cons(1); col3_uncons=col3_cons(1); %ensures only one selection occurs in case there is a resolution issue
R3uncoord_star_cons=R3_vector(col3_uncons); %optimal rotation length in community 2 under conservative non-coordination
R3_uncoord_assump_cons=R2_vector(row3_uncons); %community 2 assumption of neighbor's choice (i.e., that they will not nourish)

Neff_uncoord_cons=N_efficiency_storage(row2_uncons,col3_uncons); %system's nourishment efficiency under risky non-coordination
Neff_uncoord_risk=N_efficiency_storage(row2_unrisk,col3_unrisk); %system's nourishment efficiency under conservative non-coordination
M_uncoord_risk=TNB_storage(row2_unrisk,col3_unrisk); %system's total net benefit under risky non-coordination
M_uncoord_cons=TNB_storage(row2_uncons,col3_uncons); %system's total net benefit under conservative non-coordination
NB2_uncoord_risk=NB2_storage(row2_unrisk,col3_unrisk); %community 1's net benefit under risky non-coordination
NB2_uncoord_cons=NB2_storage(row2_uncons,col3_uncons); %community 1's net benefit under conservative non-coordination
NB3_uncoord_risk=NB3_storage(row2_unrisk,col3_unrisk); %community 2's net benefit under risky non-coordination
NB3_uncoord_cons=NB3_storage(row2_uncons,col3_uncons); %community 2's net benefit under conservative non-coordination
behavior_uncoord_cons=behavior_storage(row2_uncons,col3_uncons); %emergent behaviors for optimal rotation lengths under conservative non-coordination
behavior_uncoord_risk=behavior_storage(row2_unrisk,col3_unrisk); %emergent behaviors for optimal rotation lengths under risky non-coordination

