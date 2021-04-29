function Q_v_slice=calculate_Q_v_slice(K_v,R,r_w)
Q_v_slice_all=zeros(size(R));%the vertical flow from each slice
R_tmp=R';
idx_R=find(R_tmp>(r_w));
Q_v_slice_all(idx_R)=-pi*(R_tmp(idx_R).^2-r_w^2)*K_v;
Q_v_slice_all_shift=[Q_v_slice_all(2:end);0];
%the net vertical flow from each slice, acounting for the slice above
Q_v_slice=Q_v_slice_all-Q_v_slice_all_shift;