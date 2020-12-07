tic
gamma_vec=linspace(5,10,21); phi_vec=linspace(10,50,21);
nn=length(gamma_vec);
mm=length(phi_vec);
R2_coord=zeros(nn,mm); R3_coord=zeros(nn,mm);
R2_uncoord_cons=zeros(nn,mm); R3_uncoord_cons=zeros(nn,mm);
R2_uncoord_risk=zeros(nn,mm); R3_uncoord_risk=zeros(nn,mm);
% ratio_coord=zeros(nn,mm); ratio_uncoord_cons=zeros(nn,mm); ratio_uncoord_lib=zeros(nn,mm);
TNB_coord=zeros(nn,mm); TNB_uncoord_cons=zeros(nn,mm); TNB_uncoord_risk=zeros(nn,mm);
NB2coord=zeros(nn,mm); NB3coord=zeros(nn,mm);
NB2uncoord_cons=zeros(nn,mm); NB3uncoord_cons=zeros(nn,mm);
NB2uncoord_risk=zeros(nn,mm); NB3uncoord_risk=zeros(nn,mm);
behaviorcoord_storage=zeros(nn,mm); behavioruncoord_cons_storage=zeros(nn,mm); behavioruncoord_risk_storage=zeros(nn,mm);
TNB_diff_cons=zeros(nn,mm); TNB_diff_risk=zeros(nn,mm);
Neff_coord_storage=zeros(nn,mm); Neff_uncoord_cons_storage=zeros(nn,mm); Neff_uncoord_risk_storage=zeros(nn,mm);


%% Main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:numel(gamma_vec) 
    R2_coord_vector=zeros(1,mm); R3_coord_vector=zeros(1,mm);
    R2_uncoord_cons_vector=zeros(1,mm); R3_uncoord_cons_vector=zeros(1,mm);
    R2_uncoord_risk_vector=zeros(1,mm); R3_uncoord_risk_vector=zeros(1,mm);
%     ratio_coord_vector=zeros(1,mm); ratio_uncoord_cons_vector=zeros(1,mm); ratio_uncoord_lib_vector=zeros(1,mm);
    TNB_coord_vector=zeros(1,mm); TNB_uncoord_cons_vector=zeros(1,mm);TNB_uncoord_risk_vector=zeros(1,mm);
    NB2coord_vector=zeros(1,mm); NB3coord_vector=zeros(1,mm);
    NB2uncoord_cons_vector=zeros(1,mm); NB3uncoord_cons_vector=zeros(1,mm);
    NB2uncoord_risk_vector=zeros(1,mm); NB3uncoord_risk_vector=zeros(1,mm);
    behaviorcoord_storage_vector=zeros(1,mm); behavioruncoord_cons_storage_vector=zeros(1,mm); behavioruncoord_risk_storage_vector=zeros(1,mm);
    TNB_diff_cons_vector=zeros(1,mm); TNB_diff_risk_vector=zeros(1,mm);
    Neff_coord_vector=zeros(1,mm); Neff_uncoord_cons_vector=zeros(1,mm); Neff_uncoord_risk_vector=zeros(1,mm);
  
    for jj=1:numel(phi_vec)
        gamma=gamma_vec(ii);
        phi=phi_vec(jj);
        [R2coord_star,R3coord_star,R2uncoord_star_cons,R3uncoord_star_cons,R2uncoord_star_risk,R3uncoord_star_risk,M_coord,M_uncoord_cons,M_uncoord_risk,NB2_coord,NB3_coord,NB2_uncoord_cons,NB3_uncoord_cons,NB2_uncoord_risk,NB3_uncoord_risk,behavior_coord,behavior_uncoord_cons,behavior_uncoord_risk,Neff_coord,Neff_uncoord_cons,Neff_uncoord_risk]=maincode(gamma,phi);
        %% Storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ratio_coord_vector(jj)=Vratio_coord(1); ratio_uncoord_cons_vector(jj)=Vratio_uncoord_cons; ratio_uncoord_lib_vector(jj)=Vratio_uncoord;
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
%         VL_coord_vector(jj)=VLlost_coord(1); VL_uncoord_cons_vector(jj)=VLlost_uncoord_cons; VL_uncoord_lib_vector(jj)=VLlost_uncoord;
%         VL_equal_coord_vector(jj)=VL_equaleffort_coord(1); VL_diff_vector(jj)=VL_coord_vector(jj)-VL_equal_coord_vector(jj);
    end
%     ratio_coord(ii,:)=ratio_coord_vector; ratio_uncoord_cons(ii,:)=ratio_uncoord_cons_vector; ratio_uncoord_lib(ii,:)=ratio_uncoord_lib_vector;
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
time_elapsed=toc;
save('future_erosion_sandcost_DATA');
