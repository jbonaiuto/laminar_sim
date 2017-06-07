function plot_patch_size_luzia(subj_info, session_num, freq, snr, varargin)
% PLOT_PATCH_SIZE_LUZIA  Plot patch size results as in Troebinger 2014
%
% Use as
%   plot_patch_size_luzia(subj_info, 1, [10 30], -20)
% where the first argument is the subject info structure (from create_subjects),
% the second is the session numner, the third is the frequency range of the
% simulated data (Hz), and the fourth is the SNR (db)
% 
%   plot_patch_size_luzia(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * sources - 'both' (default) or ('pial','white') - which sources to
%    plot
%    * nsims - 60 (default) or integer - number of simulations per surface
%    * dipole_moment - 10 (default) or integer - moment of simulated dipole
%    * surf_dir - directory containing subject surfaces

% Parse inputs
defaults = struct('nsims', 60, 'dipole_moment', 10,...
    'surf_dir', 'd:\pred_coding\surf');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

methodnames={'EBB','IID','COH','MSP'}; %% just 1 method for now
methods_to_plot=[1 4];
Nmeth=length(methodnames);

white_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_white.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,...
    sprintf('%s%s-synth', subj_info.subj_id, subj_info.birth_date),'surf',...
    'ds_pial.hires.deformed.surf.gii');

allmeshes=strvcat(white_mesh, pial_mesh);
Nmesh=size(allmeshes,1);

% Sim/recon patch size combinations
sim_patch_sizes=[5 5 10 10];
reconstruct_patch_sizes=[5 10 5 10];

% Ftrue-Fother
f_trueOther=zeros(length(methodnames),length(sim_patch_sizes),params.nsims*Nmesh);

for idx=1:length(sim_patch_sizes)
    sim_patch_size=sim_patch_sizes(idx);
    reconstruct_patch_size=reconstruct_patch_sizes(idx);

    % Load whole brain - free energy results
    fname=sprintf('allcrossF_f%d_%d_SNR%d_dipolemoment%d_sim%d_reconstruct%d.mat',...
        freq(1),freq(2),snr,params.dipole_moment, sim_patch_size,...
        reconstruct_patch_size);
    data_file=fullfile('D:\layer_sim\results\',subj_info.subj_id,...
        num2str(session_num),...
        fname);
    load(data_file);
    
    for methind=1:Nmeth,       
        for simmeshind=1:Nmesh,                            
            % F reconstructed on true - reconstructed on other
            trueF=trueF;
            otherF=squeeze(allcrossF(simmeshind,1:params.nsims,2-simmeshind+1,methind));
            f_trueOther(methind,idx,(simmeshind-1)*params.nsims+1:simmeshind*params.nsims)=trueF-otherF;
        end        
    end
end

for i=1:length(methods_to_plot)
    methind=methods_to_plot(i);
    method=methodnames{methind};
    figure();
    ax1=subplot(2,2,1);
    bar([squeeze(mean(f_trueOther(methind,4,params.nsims+1:Nmesh*params.nsims),3)) squeeze(mean(f_trueOther(methind,3,params.nsims+1:Nmesh*params.nsims),3))]);
    set(gca,'XTickLabel',{'True size (10mm)','Incorrect size (5mm)'});
    ylabel('True patch size=10mm');
    title('Simulated on superficial surface');
    ax1_ylim=get(ax1,'ylim');
    
    
    ax2=subplot(2,2,2);
    bar([squeeze(mean(f_trueOther(methind,4,1:params.nsims),3)) squeeze(mean(f_trueOther(methind,3,1:params.nsims),3))]);
    set(gca,'XTickLabel',{'True size (10mm)','Incorrect size (5mm)'});
    title('Simulated on deep surface');
    ax2_ylim=get(ax2,'ylim');
    ylim=[min([ax1_ylim ax2_ylim]) max([ax1_ylim ax2_ylim])];
    set(ax1,'ylim',ylim);
    set(ax2,'ylim',ylim);
    
    ax1=subplot(2,2,3);
    bar([squeeze(mean(f_trueOther(methind,1,params.nsims+1:Nmesh*params.nsims),3)) squeeze(mean(f_trueOther(methind,2,params.nsims+1:Nmesh*params.nsims),3))]);
    set(gca,'XTickLabel',{'True size (5mm)','Incorrect size (10mm)'});
    ylabel('True patch size=5mm');
    title('Simulated on superficial surface');
    ax1_ylim=get(ax1,'ylim');
    
    ax2=subplot(2,2,4);
    bar([squeeze(mean(f_trueOther(methind,1,1:params.nsims),3)) squeeze(mean(f_trueOther(methind,2,1:params.nsims),3))]);
    set(gca,'XTickLabel',{'True size (5mm)','Incorrect size (10mm)'});
    title('Simulated on deep surface');
    ax2_ylim=get(ax2,'ylim');
    ylim=[min([ax1_ylim ax2_ylim]) max([ax1_ylim ax2_ylim])];
    set(ax1,'ylim',ylim);
    set(ax2,'ylim',ylim);
end
