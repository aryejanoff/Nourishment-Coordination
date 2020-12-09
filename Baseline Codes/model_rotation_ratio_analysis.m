%NOTE: turn on the first line below if you want to explore the rotation
%length ratios for the base case (i.e., s = 1.5 km, Figure 8a); turn on the 
%second line below if you want to explore the model rotation length ratios 
%for the larger community size case (i.e., s = 10 km, Figure 8b)
load('basePV1_basePV2_sensitivity_analyses_DATA.mat')
% load('basePV1_basePV2_sensitivity_analyses_DATA_largerS.mat')

nn=length(a2_vec);
mm=length(a3_vec);
PV_ratio=NaN(nn,mm); %property value ratios (same as beachfront wealth ratios b/c system is symmetric
Rotc_ratio=NaN(nn,mm); %coordinated model rotation length ratios from output data file
Rotnc_ratio=NaN(nn,mm); %uncoordinated model rotation length ratios from output data file

for i=1:length(a2_vec)
    for j=1:length(a3_vec)
        if a3_vec(j)>=a2_vec(i)
            PV_ratio(i,j)=(a3_vec(j)*w_init^beta)./(a2_vec(i)*w_init^beta);
            Rotc_ratio(i,j)=R3_coord(i,j)./R2_coord(i,j);
            Rotnc_ratio(i,j)=R3_uncoord_cons(i,j)./R2_uncoord_cons(i,j);
        elseif a3_vec(j)<a2_vec(i)
            PV_ratio(i,j)=(a2_vec(i)*w_init^beta)./(a3_vec(j)*w_init^beta);
            Rotc_ratio(i,j)=R2_coord(i,j)./R3_coord(i,j);
            Rotnc_ratio(i,j)=R2_uncoord_cons(i,j)./R3_uncoord_cons(i,j);
        end
    end
end

%% Figures

figure (1)
%NOTE: turn on the first two lines below if you want to explore the rotation
%length ratios for the base case (i.e., s = 1.5 km, Figure 8a); turn on the 
%third and fourth lines below if you want to explore the model rotation length ratios 
%for the larger community size case (i.e., s = 10 km, Figure 8b)
subplot(1,2,1)
title('Smaller Communities (s = 1.5 km)') 
% subplot(1,2,2)
% title('Larger Communities (s = 10 km)') 
hold on
box on
plot(PV_ratio,Rotc_ratio,'kx')
plot(PV_ratio,Rotnc_ratio,'ro') 
xlabel('Beachfront Wealth Ratio ($$/$)')
ylabel('Rotation Length Ratio ($$/$)')
xlim([1 5])
ylim([0 1])
pbaspect([1 1 1])
set(gca,'FontSize',12)

