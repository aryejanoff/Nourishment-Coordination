load('basePV1_basePV2_sensitivity_analyses_DATA.mat')

R2_cons_diff=NaN(nn,mm);
R3_cons_diff=NaN(nn,mm);
R2_risk_diff=NaN(nn,mm);
R3_risk_diff=NaN(nn,mm);
R2_cons=NaN(nn,mm);
R3_cons=NaN(nn,mm);
R2_risk=NaN(nn,mm);
R3_risk=NaN(nn,mm);

for ii=1:numel(a2_vec) 
    for jj=1:numel(a3_vec)
        if isnan(R2_coord(ii,jj))
            R2_coord(ii,jj)=100;
        end
        if isnan(R3_coord(ii,jj))
            R3_coord(ii,jj)=100;
        end
        if isnan(R2_uncoord_cons(ii,jj))
            R2_uncoord_cons(ii,jj)=100;
        end
        if isnan(R3_uncoord_cons(ii,jj))
            R3_uncoord_cons(ii,jj)=100;
        end
        if isnan(R2_uncoord_risk(ii,jj))
            R2_uncoord_risk(ii,jj)=100;
        end
        if isnan(R3_uncoord_risk(ii,jj))
            R3_uncoord_risk(ii,jj)=100;
        end
        R2_cons_diff(ii,jj)=R2_uncoord_cons(ii,jj)-R2_coord(ii,jj);
        R3_cons_diff(ii,jj)=R3_uncoord_cons(ii,jj)-R3_coord(ii,jj);
        R2_risk_diff(ii,jj)=R2_uncoord_risk(ii,jj)-R2_coord(ii,jj);
        R3_risk_diff(ii,jj)=R3_uncoord_risk(ii,jj)-R3_coord(ii,jj);
        if R2_cons_diff(ii,jj)>0 && TNB_diff_cons(ii,jj)>0
            R2_cons(ii,jj)=-1;
        elseif R2_cons_diff(ii,jj)<0 && TNB_diff_cons(ii,jj)>0
            R2_cons(ii,jj)=1;
        end
        if R3_cons_diff(ii,jj)>0 && TNB_diff_cons(ii,jj)>0
            R3_cons(ii,jj)=-1;
        elseif R3_cons_diff(ii,jj)<0 && TNB_diff_cons(ii,jj)>0
            R3_cons(ii,jj)=1;
        end
        if R2_risk_diff(ii,jj)>0 && TNB_diff_risk(ii,jj)>0
            R2_risk(ii,jj)=-1;
        elseif R2_risk_diff(ii,jj)<0 && TNB_diff_risk(ii,jj)>0
            R2_risk(ii,jj)=1;
        end
        if R3_risk_diff(ii,jj)>0 && TNB_diff_risk(ii,jj)>0
            R3_risk(ii,jj)=-1;
        elseif R3_risk_diff(ii,jj)<0 && TNB_diff_risk(ii,jj)>0
            R3_risk(ii,jj)=1;
        end
    end
end

%% Figures
figure (1) %figure 5f-g in paper
subplot(2,1,1)
pcolor(a3_vec/1e3,a2_vec/1e3,R2_cons)
xlabel('Baseline PV_2 ($10^3)')
ylabel({'Community 1';'Baseline PV_1 ($10^3)'})
title({'Uncoordinated nourishment';'effort relative to coordination'})
set(gca,'FontSize',20)
caxis([-1 1])
pbaspect([1 1 1])
shading flat
colormap(flipud(copper))
subplot(2,1,2)
pcolor(a3_vec/1e3,a2_vec/1e3,R3_cons)
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^3)')
ylabel({'Community 2';'Baseline PV_1 ($10^3)'})
set(gca,'FontSize',20)
caxis([-1 1])
colormap(flipud(copper))

figure (2) %figure A1d-e in paper
subplot(2,1,1)
pcolor(a3_vec/1e3,a2_vec/1e3,R2_risk)
xlabel('Baseline PV_2 ($10^3)')
ylabel({'Community 1';'Baseline PV_1 ($10^3)'})
title({'Uncoordinated nourishment';'effort relative to coordination'})
set(gca,'FontSize',20)
caxis([-1 1])
pbaspect([1 1 1])
shading flat
colormap(flipud(copper))
subplot(2,1,2)
pcolor(a3_vec/1e3,a2_vec/1e3,R3_risk)
shading flat
pbaspect([1 1 1])
xlabel('Baseline PV_2 ($10^3)')
ylabel({'Community 2';'Baseline PV_1 ($10^3)'})
set(gca,'FontSize',20)
caxis([-1 1])
colormap(flipud(copper))

