function [R_plot,Z_plot,zspan,hw_all,t]=recharge_well_computation_given_Qw_wEXP_wHYD_imp(K,I,L_w,T,r_w,n,theta_i,T_on,Qw_on,eps_R,eps_h,time_int,Nz,RK,PSI_i)
tic;
% parameters (units are m, minute)
% K = horizontal hydraulic coductivity [m/min]; I = anisotropy ration [-]; L_w = well's length [m]; T = total simulation time [days]
% r_w = well's radius [m]; n = porosity [-]; theta_i = initail water content [-]; T_on = time in which injection is on [min]; 
% Qw_on = discharge into the well when the discharge is on [m^3/min];
% eps_R = wetting front location tolarence [-]; eps_h = well water level
% tolarence [-]; time_int = number of time intervals [-]; RK = number of
% Runge–Kutta sections [-]

 K_h = K;          
 K_v = K/I;     
 
 zspan=linspace(0,L_w,Nz); 
 A=pi*r_w^2;
 WI = 1.01-PSI_i/(L_w-PSI_i);
 dt_max = A/(10*2*pi*K_h*(L_w-PSI_i));
 t=linspace(0,T,time_int);
 dt_all=diff(t);
 dz_max = (T_on*Qw_on)/(10*2*pi*r_w*(L_w-PSI_i));
 
 if dt_all(1)>dt_max
     time_int = ceil(T/dt_max);
     t=linspace(0,T,time_int);
    dt_all=diff(t);
 end
 
 if zspan(2)>dz_max
     Nz = ceil(L_w/dz_max);
     zspan=linspace(0,L_w,Nz); 
 end

 Qw_all=zeros(size(t));
 hw_all=zeros(size(t));
 Q_spill=zeros(size(t));
 V_total=zeros(size(t));
 V_pm=zeros(size(t)); 
 Qw_all(t<=T_on)=Qw_on;
 R_all=zeros(length(zspan),length(t)); 
 R_all(:,1)=r_w; % initial wetting front is at r=r_w, i.e. no wetting.
 Z_vertical=zeros(size(t));
 R_vertical=r_w*ones(size(t));

 for ti=1:numel(t)-1 
     dt=dt_all(ti); % present dt
     temp=Qw_all(ti)+Qw_all(ti+1);
     Qw= temp / 2; % the present Qw
     hw=hw_all(ti); % the present hw
     R=R_all(:,ti);
     Z=Z_vertical(ti);
     [next_hw,next_R,next_Q_spill,next_Z,next_zR,dt2]=compute_imp(dt,Qw,A,K_h,K_v,hw,R,r_w,zspan,L_w,n,theta_i,Z,PSI_i,WI,eps_R,eps_h,RK);
     if dt2<dt
         dt2_2 = dt2;
         while dt2_2<=dt
             [next_hw,next_R,next_Q_spill,next_Z,next_zR,dt2]=compute_imp(dt2,Qw,A,K_h,K_v,next_hw,next_R,r_w,zspan,L_w,n,theta_i,next_Z,PSI_i,WI,eps_R,eps_h,RK);
             dt2_2 = dt2_2+dt2;
         end
     end
     
     hw_all(ti+1)=next_hw;
     R_all(:,ti+1)=next_R;
     Q_spill(ti+1)=next_Q_spill; 
     Z_vertical(ti+1)=next_Z;
     R_vertical(ti+1)=next_zR;
     V_total(ti+1)=V_total(ti)+Qw*dt; % total volume of water that was recharged into the well
     V_pm(ti+1)=V_total(ti+1)-next_hw*A; % water volume in porous media is total vol. minus vol. in well.
 end
 
 Z_vert_origin=Z_vertical-max(Z_vertical);
 %construct the vertical flow as a reconstruction of the bottom R
 R_v=R_all(1,:);

% t_to_show=[2,5, 10, 20, 30, 40, 60, 80, 120, 180]; % times to show in figure

h_w=hw_all;
Z_vert_fl=flip(Z_vert_origin)';
R_vert_fl=r_w*ones(length(t),length(t_to_show));

idx_of_t_toshow=zeros(size(t_to_show));
for i=1:numel(t_to_show)
    idx_of_t_toshow(i)=find(t>=t_to_show(i),1);
    R_vert_fl_tmp=flip(R_vertical(1:idx_of_t_toshow(i)));
    R_vert_fl(1:idx_of_t_toshow(i),i)=R_vert_fl_tmp';
end

% figure(1)
% plot(t,hw_all)

% figure(2)
 R_plot=[R_vert_fl;R_all(:,idx_of_t_toshow)];
 R_plot(R_plot==r_w)=nan;
 Z_plot=[Z_vert_fl',zspan];
%  plot(R_plot,Z_plot)
 