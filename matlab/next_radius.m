function [O_r, next_R] = next_radius(K_h, K_v, next_hw, hw, zspan, next_R, R, r_w, r_initial,PSI, n, theta_i, dt, temp_R, O_r, eps_R, RG); 

    [Q_h_slice,next_R]=calculate_Q_slice(K_h,(next_hw+hw)/2,zspan,(next_R+R)./2,r_w,r_initial,PSI); 
    Q_v_slice=calculate_Q_v_slice(K_v,(next_R+R)./2,r_w)';
    
    next_R=R+(Q_h_slice+Q_v_slice./zspan(2))'./(2*pi*((next_R+R)./2)*(n-theta_i))*dt;
    
    O_r(2) = mean(abs(next_R-temp_R)./R);
    
    if O_r(2)<O_r(1)
        O_r(1) = O_r(2);
    else
        while  O_r(2)>eps_R & O_r(2)>O_r(1)
            temp_R = next_R;
            
            next_R_Q=R+(Q_h_slice+Q_v_slice./zspan(2))'./(2*pi*(R)*(n-theta_i))*(dt/RG);
            R_T = R;
            
            for i = 2:RG-1
            next_R_Q=R_T+(Q_h_slice+Q_v_slice./zspan(2))'./(2*pi*((R_T+next_R_Q)./2)*(n-theta_i))*(dt/RG);
            R_T = next_R_Q;
            end
            
            next_R=R_T+(Q_h_slice+Q_v_slice./zspan(2))'./(2*pi*((next_R+R_T)./2)*(n-theta_i))*(dt/RG);
            
            O_r(2) = mean(abs(next_R-temp_R)./R);
        end
            [Q_h_slice,next_R]=calculate_Q_slice(K_h,(next_hw+hw)/2,zspan,(next_R+R)./2,r_w,r_initial,PSI); 
            Q_v_slice=calculate_Q_v_slice(K_v,(next_R+R)./2,r_w)';
            
            next_R=R+(Q_h_slice-Q_v_slice./zspan(2))'./(2*pi*((next_R+R)./2)*(n-theta_i))*dt;
            
            O_r(1) = mean(abs(next_R-temp_R)./R);
    end
end