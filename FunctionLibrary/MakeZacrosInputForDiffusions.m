function MakeZacrosInputForDiffusions(PES,MinsLoc,Barriers)

if nargin==1
    MinsLoc=PES.Mins;
    Barriers=PES.Barriers;
    PES=PES.Class;
end

PES=PES.interpolatePESGrid;

sitenames=1:length(MinsLoc);
NanList=cellfun(@(V) any(isnan(V(:))), Barriers.Grid);
EmptyList=cellfun(@(V) any(isempty(V(:))), Barriers.Grid);
NeighborList=1-NanList-EmptyList;
NeighborMat=(NeighborList.*(sitenames'*ones(1,length(MinsLoc))))';

Lattice.Name=sitenames;
Lattice.Position=MinsLoc;
Lattice.Type=sitenames;

Lattice.Neighbors=mat2cell(NeighborMat,[ones(1,length(MinsLoc))]);
maxcoordination=1;
for i=1:length(MinsLoc)
    Lattice.Neighbors{i}=Lattice.Neighbors{i}(Lattice.Neighbors{i}~=0);
    Lattice.Coordination{i}=length(Lattice.Neighbors{i});
    Lattice.Energy{i}=getMEPEnergies(PES,MinsLoc(i,1),MinsLoc(i,2));
    maxcoordination=max([maxcoordination Lattice.Coordination{i}]);
end

fid = fopen ('ZacrosInputs/energetics_input.dat', 'wt' );
fprintf(fid, '# Energetics for H diffusion \n\n energetics \n\n ############################################################################ \n\n');
for i=1:length(MinsLoc)
    fprintf(fid, 'cluster H \n   sites 1 \n   lattice_state \n   1 H*  1 \n\n');
    fprintf(fid, 'variant %.0f # ClusterNum%.0f \n',Lattice.Name(i),i);
    fprintf(fid, '  site_types %.0f \n', Lattice.Type(i));
    fprintf(fid, '  cluster_eng %f \n', cell2mat(Lattice.Energy(i)));
    fprintf(fid, '  graph_multiplicity 1 \n');
    fprintf(fid, 'end_variant \n\n');
    fprintf(fid, 'end_cluster \n##########\n ');
end
fclose(fid);

fid = fopen ('ZacrosInputs/lattice_input.dat', 'wt' );
fprintf(fid, 'lattice explicit \n\n');
fprintf(fid, 'n_sites %.0f \n',length(MinsLoc));
fprintf(fid, 'max_coord %.0f \n',maxcoordination);
fprintf(fid, 'n_site_types %.0f \n',length(MinsLoc));
sitenamesstr=mat2str(sitenames);
fprintf(fid, 'site_type_names %s \n\n',sitenamesstr([2:end-1]));
fprintf(fid, 'lattice_structure \n');
for i=1:length(MinsLoc)
    NNstr=mat2str(Lattice.Neighbors{i});
    if length(Lattice.Neighbors{i})>1 %Need this bit of logic because mat2str treats mats different than numbers
        fprintf(fid, '%.0f %f %f %.0f %.0f %s \n',Lattice.Name(i),MinsLoc(i,1),MinsLoc(i,2),Lattice.Type(i),Lattice.Coordination{i},NNstr([2:end-1]));
    elseif length(Lattice.Neighbors{i})==1
        fprintf(fid, '%.0f %f %f %.0f %.0f %s \n',Lattice.Name(i),MinsLoc(i,1),MinsLoc(i,2),Lattice.Type(i),Lattice.Coordination{i},NNstr);
    end
end
fprintf(fid, 'end_lattice_structure \n');
fprintf(fid, 'end_lattice \n');
fclose(fid);

fid = fopen ('ZacrosInputs/mechanism_input.dat', 'wt' );
fprintf(fid, '# HCOOH diffusion on Au Cluster \n\n mechanism \n\n ############################################################################ \n\n');

for i=1:length(Barriers.Name)
    
    %Check that it is the right value in the energies
    if ~(cell2mat(Lattice.Energy(Barriers.StartMinIndex{i}))==Barriers.Energies{i}(1))||...
            ~(cell2mat(Lattice.Energy(Barriers.EndMinIndex{i}))==Barriers.Energies{i}(3))
        error('Error \nThe MEP minimums for %d does not match the minimums',i)
    end
    
    fprintf(fid, 'reversible_step diff/H/%s \n',Barriers.Name{i});
    fprintf(fid, '  sites 2 \n');
    fprintf(fid, '  neighboring 1-2 \n');
    fprintf(fid, '  initial \n');
    fprintf(fid, '  1 H*  1 \n');
    fprintf(fid, '  2  *  1 \n');
    fprintf(fid, '  final \n');
    fprintf(fid, '  2  *  1 \n');
    fprintf(fid, '  1 H*  1 \n\n');
    
    fprintf(fid, 'variant %s VariantNum1\n',Barriers.Name{i});
    fprintf(fid, '  site_types %.0f %.0f \n',Barriers.StartMinIndex{i},Barriers.EndMinIndex{i});
    fprintf(fid, '  pre_expon 1.00e+13 \n');
    fprintf(fid, '  pe_ratio 1.00e+00 \n');
    fprintf(fid, '  activ_eng %f \n',Barriers.Energies{i}(2)-Barriers.Energies{i}(1));
    fprintf(fid, '  prox_factor 0.50 \n');
    fprintf(fid, '#  stiffness_scalable \n');
    fprintf(fid, 'end_variant \n\n');
    fprintf(fid, 'end_reversible_step \n########## \n');
end
fprintf(fid, 'end_mechanism');
fclose(fid);

fid = fopen ('ZacrosInputs/simulation_input.dat', 'wt' );
fprintf(fid, '# KMC simulation specification \n\n');
fprintf(fid, 'random_seed               123 \n\n');
fprintf(fid, 'temperature               413.00 \n');
fprintf(fid, 'pressure                  1.00 \n\n');
fprintf(fid, 'n_gas_species             1 \n');
fprintf(fid, 'gas_specs_names           H2 \n');
fprintf(fid, 'gas_energies              0.000 \n');
fprintf(fid, 'gas_molec_weights         2.016 \n');
fprintf(fid, 'gas_molar_fracs           1.000 \n');
fprintf(fid, 'n_surf_species            1 \n');
fprintf(fid, 'surf_specs_names          H* \n');
fprintf(fid, 'surf_specs_dent           1 \n\n');
fprintf(fid, 'snapshots                 on event 100000000 \n');
fprintf(fid, 'process_statistics        on event 1000000 \n');
fprintf(fid, 'species_numbers           on event 100000000 \n\n');
fprintf(fid, 'event_report              on elemevent 1000000 \n\n');
fprintf(fid, 'max_steps                 100000000 \n');
fprintf(fid, 'max_time                  100 \n\n');
fprintf(fid, 'wall_time                 43200 \n\n');
fprintf(fid, '# no_restart \n\n');

fprintf(fid, '# enable_stiffness_scaling \n');
fprintf(fid, '# check_every 200 \n');
fprintf(fid, '# stiffness_scale_all \n');
fprintf(fid, '# min_separation 40.0 \n');
fprintf(fid, '# max_separation 60.0 \n');
fprintf(fid, '# max_qequil_separation 3.0 \n');
fprintf(fid, '# stiffn_coeff_threshold 0.10 \n');
fprintf(fid, '# energetics_list         on event \n');
fprintf(fid, '# debug_report_processes \n');
fprintf(fid, '# debug_report_global_energetics \n');
fprintf(fid, '# debug_check_processes \n');
fprintf(fid, '# debug_check_lattice \n\n');
fprintf(fid, 'finish \n\n');
fclose(fid);

fid = fopen ('ZacrosInputs/state_input.dat', 'wt' );
fprintf(fid, 'initial_state \n');
fprintf(fid, 'seed_on_sites H* 1 \n');
fprintf(fid, 'end_initial_state \n');
fclose(fid);
end