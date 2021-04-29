function [next_hw,next_R,next_Q_spill,next_Z,next_zR,dt]=compute_imp(dt,Qw,A,K_h,K_v,hw,R,r_w,zspan,L_w,n,theta_i,Z,PSI,WI,eps_R,eps_h,RK)
% compute water level in well, position of wetting front, and spilled
% volume at a time step
% as  a function of the present variables: 
%   - R(z) the position of wetting front, 
%   - Qw the discharge into the well
%   - dt the present tiem step
%
% first, calculate Q_slice(z) for all z values
N = 0;
O_r = 1;
O_h = 1;
r_initial = WI;
next_R = R;
next_hw = hw;
while O_r(1)>eps_R | O_h(1)>eps_h
    temp_R = next_R;
    temp_hw = next_hw;
        
    [O_r, next_R] = next_radius(K_h,K_v,next_hw,hw,zspan,next_R,R,r_w,r_initial,PSI,n,theta_i,dt,temp_R,O_r,eps_R,RK); 
       
    [O_h, next_hw,  next_Q_spill, next_Z, next_zR] = next_hight(K_v,n,theta_i,dt,Z,next_R,R,zspan,r_w,hw,Qw,A,L_w,temp_hw,O_h); 


    if O_h(2)>=O_h(1) || next_hw < 0
        dt = dt/2;
        O_h(1) = 1;
    else

        O_h(1) = O_h(2);
    end
end