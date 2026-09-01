% This code is used to visualize ModEM models
% planes can be defined which are plotted over the model

m = load_model_modem('DIKE1_NLCG_021.rho'); %loads model
dobs = load_data_modem('Datafile.data'); %loads observed data
[m,dobs] = link_model_data(m,dobs); %geo-references model
%plot_diagonal_section(m,dobs); %plot diagonal slice

% planes in format [NS1 EW1 NS2 EW2,...]
planes = [0.05 -0.013 -0.049 0.0033; -0.0075 0.0198 -0.011 -0.048];
plot_model_3D_RV_planes(m,dobs,planes)
plot_model_3D(m,dobs)

%plot_diagonal_section(m,dobs)

plot_slice(m,15,dobs)