function [O_h, next_hw,  next_Q_spill, next_Z, next_zR] = next_hight(K_v,n,theta_i,dt,Z,next_R,R,zspan,r_w,hw,Qw,A,L_w,temp_hw,O_h)
     

    dz=K_v/(n-theta_i)*dt;
    
    next_Z=Z+dz;
    next_zR=next_R(1);

    Q_pm = sum(pi.*(next_R.^2-R.^2).*(n-theta_i).*zspan(2))+dz*pi*(next_zR^2-r_w^2)*(n-theta_i);
    if Q_pm<0
        Q_pm=0;
    end
    next_hw_theoretical=hw+(Qw*dt-Q_pm)/A;
    
    if next_hw_theoretical>L_w
        next_Q_spill=(next_hw_theoretical-L_w)*A/dt; 
        next_hw=L_w; 
    else
        next_Q_spill=0;
        next_hw=next_hw_theoretical;
        if next_hw<0
            next_hw = 0;
        end
    end
    
        O_h(2) = abs(next_hw-temp_hw)/L_w;
            
end
