
a=H_SingleSlit(512,512/2+1);

function H_Single=H_SingleSlit (z,dim,center)
H_Single = zeros(dim,dim);
v_single_slit_width   = 100;   % width pixels
H_SingleSlit_hight  = 10;    %  hight pixels

H_Single((center-H_SingleSlit_hight):(center+H_SingleSlit_hight),(center-v_single_slit_width):(center+v_single_slit_width)) =ones(2*H_SingleSlit_hight+1,2*v_single_slit_width+1);
end

 


