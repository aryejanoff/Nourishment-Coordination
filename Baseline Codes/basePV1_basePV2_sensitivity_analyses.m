tic
%% Input parameters %%
w_init=30; %initial beach width (m)
beta=0.4; %hedonic beach parameter
PV_min=0.1e6; %minimum property value ($)
PV_max=2e6; %maximum property value ($)
PV=linspace(PV_min,PV_max,20); %property value vector
a2_vec=PV/(w_init^beta); a3_vec=PV/(w_init^beta); %baseline property value vector for communities 1 and 2
nn=length(a2_vec); mm=length(a3_vec); %vector lengths

%% variable storage %%
R2_coord=zeros(nn,mm); R3_coord=zeros(nn,mm); %coordinated rotation lengths (yrs)
R2_uncoord_cons=zeros(nn,mm); R3_uncoord_cons=zeros(nn,mm); %conservative uncoordinated rotation lengths (yrs)
R2_uncoord_risk=zeros(nn,mm); R3_uncoord_risk=zeros(nn,mm); %risky uncoordinated rotation lengths (yrs)
TNB_coord=zeros(nn,mm); TNB_uncoord_cons=zeros(nn,mm); TNB_uncoord_risk=zeros(nn,mm); %total net benefits for coordination; conservative non-coordination; risky non-coordination ($)
NB2coord=zeros(nn,mm); NB3coord=zeros(nn,mm); %net benefits for communities 1 and 2 under coordination ($)
NB2uncoord_cons=zeros(nn,mm); NB3uncoord_cons=zeros(nn,mm); %net benefits for communities 1 and 2 under conservative non-coordination ($)
NB2uncoord_risk=zeros(nn,mm); NB3uncoord_risk=zeros(nn,mm); %net benefits for communities 1 and 2 under risky non-coordination ($)
behaviorcoord_storage=zeros(nn,mm); behavioruncoord_cons_storage=zeros(nn,mm); behavioruncoord_risk_storage=zeros(nn,mm); %emergent system behaviors under coordination; conservative non-coordination; risky non-coordination
TNB_diff_cons=zeros(nn,mm); TNB_diff_risk=zeros(nn,mm); %difference b/w coordinated and uncoordinated total net benefits (i.e., coordination incentive) under conservative and risky non-coordination ($)
Neff_coord_storage=zeros(nn,mm); Neff_uncoord_cons_storage=zeros(nn,mm); Neff_uncoord_risk_storage=zeros(nn,mm); %physical nourishment efficiencies for coordination; conservative non-coordination; risky non-coordination

%% Main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:numel(a2_vec) %baseline PV1 for loop
    R2_coord_vector=zeros(1,mm); R3_coord_vector=zeros(1,mm);
    R2_uncoord_cons_vector=zeros(1,mm); R3_uncoord_cons_vector=zeros(1,mm);
    R2_uncoord_risk_vector=zeros(1,mm); R3_uncoord_risk_vector=zeros(1,mm);
    TNB_coord_vector=zeros(1,mm); TNB_uncoord_cons_vector=zeros(1,mm);TNB_uncoord_risk_vector=zeros(1,mm);
    NB2coord_vector=zeros(1,mm); NB3coord_vector=zeros(1,mm);
    NB2uncoord_cons_vector=zeros(1,mm); NB3uncoord_cons_vector=zeros(1,mm);
    NB2uncoord_risk_vector=zeros(1,mm); NB3uncoord_risk_vector=zeros(1,mm);
    behaviorcoord_storage_vector=zeros(1,mm); behavioruncoord_cons_storage_vector=zeros(1,mm); behavioruncoord_risk_storage_vector=zeros(1,mm);
    TNB_diff_cons_vector=zeros(1,mm); TNB_diff_risk_vector=zeros(1,mm);
    Neff_coord_vector=zeros(1,mm); Neff_uncoord_cons_vector=zeros(1,mm); Neff_uncoord_risk_vector=zeros(1,mm);
  
    for jj=1:numel(a3_vec) %baseline PV2 for loop
        alpha2=a2_vec(ii); %base PV1 value for each element in vector
        alpha3=a3_vec(jj); %base PV2 value for each element in vector
        [R2coord_star,R3coord_star,R2uncoord_star_cons,R3uncoord_star_cons,R2uncoord_star_risk,R3uncoord_star_risk,M_coord,M_uncoord_cons,M_uncoord_risk,NB2_coord,NB3_coord,NB2_uncoord_cons,NB3_uncoord_cons,NB2_uncoord_risk,NB3_uncoord_risk,behavior_coord,behavior_uncoord_cons,behavior_uncoord_risk,Neff_coord,Neff_uncoord_cons,Neff_uncoord_risk]=angle_maincode(alpha2,alpha3); %[what we want from the maincode]=maincode(what we feed the maincode);
        %% Storage (see variable storage section for descriptions) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R2_coord_vector(jj)=R2coord_star;
        R3_coord_vector(jj)=R3coord_star;
        R2_uncoord_cons_vector(jj)=R2uncoord_star_cons;
        R3_uncoord_cons_vector(jj)=R3uncoord_star_cons;
        R2_uncoord_risk_vector(jj)=R2uncoord_star_risk;
        R3_uncoord_risk_vector(jj)=R3uncoord_star_risk;
        TNB_coord_vector(jj)=M_coord;
        TNB_uncoord_cons_vector(jj)=M_uncoord_cons;
        TNB_uncoord_risk_vector(jj)=M_uncoord_risk;
        NB2coord_vector(jj)=NB2_coord;
        NB3coord_vector(jj)=NB3_coord;
        NB2uncoord_cons_vector(jj)=NB2_uncoord_cons;
        NB3uncoord_cons_vector(jj)=NB3_uncoord_cons;
        NB2uncoord_risk_vector(jj)=NB2_uncoord_risk;
        NB3uncoord_risk_vector(jj)=NB3_uncoord_risk;
        behaviorcoord_storage_vector(jj)=behavior_coord;
        behavioruncoord_cons_storage_vector(jj)=behavior_uncoord_cons;
        behavioruncoord_risk_storage_vector(jj)=behavior_uncoord_risk;
        Neff_coord_vector(jj)=Neff_coord;
        Neff_uncoord_cons_vector(jj)=Neff_uncoord_cons;
        Neff_uncoord_risk_vector(jj)=Neff_uncoord_risk;
        TNB_diff_cons_vector(jj)=TNB_coord_vector(jj)-TNB_uncoord_cons_vector(jj);        
        TNB_diff_risk_vector(jj)=TNB_coord_vector(jj)-TNB_uncoord_risk_vector(jj);
    end
    R2_coord(ii,:)=R2_coord_vector;
    R3_coord(ii,:)=R3_coord_vector;
    R2_uncoord_cons(ii,:)=R2_uncoord_cons_vector;
    R3_uncoord_cons(ii,:)=R3_uncoord_cons_vector;
    R2_uncoord_risk(ii,:)=R2_uncoord_risk_vector;
    R3_uncoord_risk(ii,:)=R3_uncoord_risk_vector;
    TNB_coord(ii,:)=TNB_coord_vector;
    TNB_uncoord_cons(ii,:)=TNB_uncoord_cons_vector;
    TNB_uncoord_risk(ii,:)=TNB_uncoord_risk_vector;
    NB2coord(ii,:)=NB2coord_vector;
    NB3coord(ii,:)=NB3coord_vector;
    NB2uncoord_cons(ii,:)=NB2uncoord_cons_vector;
    NB3uncoord_cons(ii,:)=NB3uncoord_cons_vector;
    NB2uncoord_risk(ii,:)=NB2uncoord_risk_vector;
    NB3uncoord_risk(ii,:)=NB3uncoord_risk_vector;
    behaviorcoord_storage(ii,:)=behaviorcoord_storage_vector;
    behavioruncoord_cons_storage(ii,:)=behavioruncoord_cons_storage_vector;
    behavioruncoord_risk_storage(ii,:)=behavioruncoord_risk_storage_vector;
    Neff_coord_storage(ii,:)=Neff_coord_vector;
    Neff_uncoord_cons_storage(ii,:)=Neff_uncoord_cons_vector;
    Neff_uncoord_risk_storage(ii,:)=Neff_uncoord_risk_vector;
    TNB_diff_cons(ii,:)=TNB_diff_cons_vector;
    TNB_diff_risk(ii,:)=TNB_diff_risk_vector;
end

%% Save Data
time_elapsed=toc; %elapsed runtime for external script for computational accounting
save('basePV1_basePV2_sensitivity_analyses_DATA'); %saves output data (from this external script only!) in your root directory