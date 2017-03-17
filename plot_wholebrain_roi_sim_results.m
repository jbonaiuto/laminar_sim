function plot_wholebrain_roi_sim_results(subj_info, session_num, freq, snr, varargin)

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10, 'surf_dir', 'd:\pred_coding\surf',...
    'mri_dir', 'd:\pred_coding\mri', 'sim_patch_size',0, 'reconstruct_patch_size',0);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
Nmeth=length(methodnames);
allmeshes=strvcat(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii'),...
    fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii'));
Nmesh=size(allmeshes,1);
orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');

wholeBrainPialWhiteF=zeros(Nmeth,Nmesh,params.nsims);
for methind=1:Nmeth,       
    for simmeshind=1:Nmesh,    
        if params.sim_patch_size==0 && params.reconstruct_patch_size==0
            data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d.mat',freq(1),freq(2),snr,params.dipole_moment));
        else
            data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id, num2str(session_num), sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',freq(1),freq(2),snr,params.dipole_moment, params.sim_patch_size, params.reconstruct_patch_size));
        end
        load(data_file);

        % F reconstructed on true - reconstructed on other
        % num simulations x number of folds
        pialWhiteF=squeeze(allcrossF(simmeshind,1:params.nsims,2,methind)-allcrossF(simmeshind,1:params.nsims,1,methind));
        wholeBrainPialWhiteF(methind,simmeshind,:)=pialWhiteF;
    end
end

roiWhitePialT=zeros(Nmeth,Nmesh*params.nsims);
data_dir=fullfile('D:\layer_sim\ttest_results', subj_info.subj_id, num2str(session_num), sprintf('f%d_%d_SNR%d_dipolemoment%d', freq(1), freq(2), snr, params.dipole_moment));
for methind=1:Nmeth,    
    method=methodnames{methind};
    disp(method);
    
    roiWhitePialT(methind,:)=get_wmpial_t(data_dir, method, params.nsims, allmeshes(2,:), ...
        allmeshes(1,:), orig_pial_mesh, orig_white_mesh);
end

figure('Position',[1 1 1700 1800]);
mesh_idx=zeros(Nmeth,Nmesh,2);
mesh_idx(1,1,:)=[6 7];
mesh_idx(1,2,:)=[1 2];
mesh_idx(2,1,:)=[16 17];
mesh_idx(2,2,:)=[11 12];
mesh_idx(3,1,:)=[26 27];
mesh_idx(3,2,:)=[21 22];
mesh_idx(4,1,:)=[36 37];
mesh_idx(4,2,:)=[31 32];
for methind=1:Nmeth,       
    for simmeshind=1:Nmesh,    
        subplot(Nmeth*Nmesh,5,mesh_idx(methind,simmeshind,:));
        [path,file,ext]=fileparts(deblank(allmeshes(simmeshind,:)));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        hold on
        for simind=1:params.nsims
            color='b';
            if wholeBrainPialWhiteF(methind,simmeshind,simind)<0
                color='r';
            end
            bar(simind,wholeBrainPialWhiteF(methind,simmeshind,simind),color,'EdgeColor','none');
        end
        plot([0 params.nsims+1],[3 3],'k--');
        plot([0 params.nsims+1],[-3 -3],'k--');
        xlim([0 params.nsims+1]);
        xlabel('Simulation')
        ylabel('Free energy diff (pial-white)');
        title(sprintf('Free energy, %s, %s',methodnames{methind},simmeshname));        
    end

end

mesh_idx=zeros(Nmeth,Nmesh,2);
mesh_idx(1,1,:)=[8 9];
mesh_idx(1,2,:)=[3 4];
mesh_idx(2,1,:)=[18 19];
mesh_idx(2,2,:)=[13 14];
mesh_idx(3,1,:)=[28 29];
mesh_idx(3,2,:)=[23 24];
mesh_idx(4,1,:)=[38 39];
mesh_idx(4,2,:)=[33 34];
for methind=1:Nmeth,    
    method=methodnames{methind};
    disp(method);
    
    % For each simulated mesh
    for simmeshind=1:Nmesh,
        [path,file,ext]=fileparts(allmeshes(simmeshind,:));
        x=strsplit(file,'.');
        y=strsplit(x{1},'_');
        simmeshname=y{2};
        
        subplot(Nmeth*Nmesh,5,mesh_idx(methind,simmeshind,:));
        hold on
        for simind=1:params.nsims
            color='b';
            if roiWhitePialT(methind,(simmeshind-1)*params.nsims+simind)<0
                color='r';
            end
            bar(simind,roiWhitePialT(methind,(simmeshind-1)*params.nsims+simind),color,'EdgeColor','none');
        end
        hold on
        plot([0 params.nsims+1], [1.69 1.69],'k--');
        plot([0 params.nsims+1], [-1.69 -1.69],'k--');
        xlim([0 params.nsims+1]);
        %ylim([-40 40]);
        xlabel('Simulation')
        ylabel('Mean t (pial-white)');
        title(sprintf('t Stat, %s, %s',methodnames{methind},simmeshname)); 
    end
end

meth_idx=[5 10;15 20;25 30; 35 40];
for methind=1:Nmeth
    subplot(Nmeth*Nmesh,5,meth_idx(methind,:));
    hold on;
    plot(squeeze(wholeBrainPialWhiteF(methind,1,:)),roiWhitePialT(methind,1:params.nsims),'o','MarkerEdgeColor','none','MarkerFaceColor','r');
    plot(squeeze(wholeBrainPialWhiteF(methind,2,:)),roiWhitePialT(methind,params.nsims+1:end),'o','MarkerEdgeColor','none','MarkerFaceColor','b');

    f=[squeeze(wholeBrainPialWhiteF(methind,1,:)); squeeze(wholeBrainPialWhiteF(methind,2,:))];
    [rho,pval] = corr(f, roiWhitePialT(methind,:)','type','Spearman');    
    pPoly = polyfit(f, roiWhitePialT(methind,:)', 1); % Linear fit of xdata vs ydata
    linePointsX = [min(f) max(f)]; % find left and right x values
    linePointsY = polyval(pPoly,[min(f),max(f)]); 
    plot (linePointsX,linePointsY,'--k')
    hold off;
    xlabel('F');
    ylabel('t');
    title(sprintf('%s, rho=%.2f, p=%0.5f\n', methodnames{methind},rho,pval));
    legend({'Deep source','Superficial source'});
end
