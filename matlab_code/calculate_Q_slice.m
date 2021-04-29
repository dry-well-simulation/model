function [Q_slice,R_correct]=calculate_Q_slice(K,hw,zspan,R,r_w,r_initial,PSI)
R_correct = R;
Q_slice=zeros(size(zspan));
idx=find(zspan<hw);
%
if ~isempty(idx)
    if length(idx)>1
        idx=idx(1:end-1);
    else
        idx=[];
    end
end
%
for i=idx
    z=zspan(i);
    eps_x = R(i)/r_w;
    if eps_x>r_initial
        Q_slice(i)=2*pi*K*(hw-z-PSI) / log (eps_x) ;
    else
        
        eps_x = r_initial;
        R_correct(i) = r_w*eps_x;
        Q_slice(i)=2*pi*K*(hw-z-PSI) / log (eps_x) ;
    end

    if Q_slice(i)<0
        Q_slice(i)=0;
    end
end