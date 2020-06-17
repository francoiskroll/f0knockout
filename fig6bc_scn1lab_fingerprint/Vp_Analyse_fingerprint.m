%% Vp_Analyse_MD

% marcus.ghosh.11@ucl.ac.uk

%% Load

%Load Data - using multiselect
[filename, pathname] = uigetfile('*.mat', 'Select file or files','MultiSelect','on'); %Select files
if isequal(filename,0) %If no file is selected
    error('No Files Selected') %Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) %Show selected filenames
end

%% Options 

ct{1} = [2 3]; % days of interest 
ct{2} = [2 3]; % nights of interest  

%% Combining Experiments (< Parameter Extracted Data)

% Convert filename to a cell (useful for single files)
filename = cellstr(filename);

experiment_tags = [];
group_tags = [];
parameter_matrix = [];
for f = 1:size(filename,2) %For each file
    clear experiment;
    experiment = load(strcat(pathname,filename{f})); %Load the mat file
    experiment_tags = [experiment_tags ; ones(size(experiment.group_tags,1),1)*f];
    group_tags = [group_tags ; experiment.group_tags];
    % Allocate experiment tags
    
    % Load Data 
    parameter_matrix = [parameter_matrix ; experiment.parameter_matrix];
    
    % Nab variables
    if f == 1 % For the first file
        parameters = experiment.parameters;
        time_window = experiment.time_window;
        time_bins = experiment.time_bins;
        unit_conversion = experiment.unit_conversion;
        cmap = experiment.cmap;
        cmap_2 = experiment.cmap_2;
        night_color = experiment.night_color;
        nights = experiment.nights;
        nights_crop = experiment.nights_crop;
        geno_list = experiment.geno_list;
        units = experiment.units;
        units_2 = experiment.units_2;
        parameter_smooth = experiment.parameter_smooth;
        first_night = experiment.first_night;
        days_crop = experiment.days_crop;
        days = experiment.days;
    end
    
    
end

clear pathname filename geno f;

%%  Adjust Colors 

cmap_2 = flip(cmap_2); % flip cmap
for c = 1:2:size(cmap_2,1)
    cmap_2([c c+1],:) = cmap_2([c+1 c],:); % swap the colors around
end
cmap = cmap_2(1:2:size(cmap_2,1),:); % Extract main colors

%% Correct Group Tags for Francois 

% ! only if loaded in that specific order: 
% 140218_34_DATA.mat;
% 190813_06_DATA.mat;
% 190813_07_DATA.mat

scrap = group_tags(experiment_tags >= 2); 
scrap(scrap == 1) = 3; 
scrap(scrap == 2) = 1; 

group_tags(experiment_tags >= 2) = scrap; 

%% Format Data (fish x parameters)

days = ct{1}; % days of interest 
nights = ct{2}; % nights of interest 
time_window(1) = min([days_crop(days) nights_crop(nights)]);
time_window(2) = max([days_crop(days) nights_crop(nights)]);

% Parameter_matrix
scrap = nanmean(parameter_matrix(:,1:10,days_crop(days)),3); % day average % removes the last two parameters
scrap = [scrap nanmean(parameter_matrix(:,1:10,nights_crop(nights)),3)]; % combine with night average

parameter_matrix = scrap; 

%% Z-Score per experiment

for e = 1:max(experiment_tags) % for each experiment 
    wt_mean = nanmean(parameter_matrix(experiment_tags == e & group_tags == 1,:)); % wt - mean
    wt_std = nanstd(parameter_matrix(experiment_tags == e & group_tags == 1,:)); % wt - std
    
    parameter_matrix(experiment_tags == e,:) = ...
        (parameter_matrix(experiment_tags == e,:) - wt_mean)./wt_std; % z-score to wt mean
end 

%% FingerPrints 

figure; hold on;
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font
plot([0 size(parameter_matrix,2)], [0 0],'--k','linewidth',3); % baseline 
        
for e = 1:max(experiment_tags) % for each experiment
    
    for g = 3:max(group_tags(experiment_tags == e)) % plot only Hom (not Wt or Het)
        
        if e == 1 % for Ellen's experiment; in black
            errorbar(nanmean(parameter_matrix(experiment_tags == e & group_tags == g,:)),...
                nanstd(parameter_matrix(experiment_tags == e & group_tags == g,:)),...
                'color','k','linewidth',3);
        else
            errorbar(nanmean(parameter_matrix(experiment_tags == e & group_tags == g,:)),...
                nanstd(parameter_matrix(experiment_tags == e & group_tags == g,:)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);
        end
        
    end
    
end

% Night Shading 
x_lim = get(gca,'XLim');  %# Get the range of the x axis
y_lim = get(gca,'YLim');  %# Get the range of the y axis
r(1) = rectangle('Position',[10.5 y_lim(1) 9.5  sum(abs(y_lim))],...
    'FaceColor',night_color,'Edgecolor',[1 1 1]);
uistack(r(1),'bottom'); % Send to back

% Axis Labels 
ylabel('Deviation from Wild-Type'); 
xlabel('Parameters');

parametermatrix_fingerprint = parameter_matrix; %FK for export

%% Z-Score matrix 
parameter_matrix = zscore(parameter_matrix); % z-score to standardise parameters

%% Euclidean Distance - From Ellen's Data 
ED(:,1) = pdist2(parameter_matrix, ...
    nanmean(parameter_matrix(experiment_tags == 1 & group_tags == 1,:)));
ED(:,2) = pdist2(parameter_matrix, ...
    nanmean(parameter_matrix(experiment_tags == 1 & group_tags == 3,:)));

%% Scatter Plot 

figure; hold on;
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font
plot([min(ED(:,2)), max(ED(:,1))],[min(ED(:,1)), max(ED(:,1))],'--k','linewidth',3); % baseline 

for e = 1:max(experiment_tags)
    
    for g = 1:max(group_tags(experiment_tags == e))%% Distance Figure
        
        scatter(ED(experiment_tags == e & group_tags == g,2),...
            ED(experiment_tags == e & group_tags == g,1),90,...
            'markerfacecolor',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),...
            'markeredgecolor',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),...
            'markerfacealpha',0.5)

    end 
end 

axis([min(ED(:,2)), max(ED(:,1)), min(ED(:,1)), max(ED(:,1))]); 

xlabel('Distance from -/-')
ylabel('Distance from WT')

parametermatrix_scatterHom = parameter_matrix;

%% Euclidean Distance from Paired WT's 

figure;
hold on; box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); set(gca,'FontName','Calibri'); % Set Font

g_counter = 1;
for e = 1:max(experiment_tags)
    for g = 1:max(group_tags)
        spread_cols = plotSpread(...
            pdist2(parameter_matrix(experiment_tags == e & group_tags == g,:),...
            nanmean(parameter_matrix(experiment_tags == e & group_tags == 1,:))),...
            'xValues',g_counter,...
            'distributionColors',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'showMM',2);
        spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
        if sum(experiment_tags == e & group_tags == g) > 0
            g_counter = g_counter + 1;
        end
    end
end

set(findall(gca,'type','line'),'markersize',30); % change

% Nice Figure
ylabel('Euclidean Distance From WT Mean');
%set(gca,'xtick',1:max(groups) + 2); % Set x ticks
%set(gca,'xticklabel',['WT Day','WT Night',geno_list{6}.colheaders]); % Name each group


% FK below is to write euclidian distances matrix for exporting
eucdists = NaN([48 9]);

counter = 1
for e = 1:max(experiment_tags)
    for g = 1:max(group_tags)
        
        scrap = pdist2(parameter_matrix(experiment_tags == e & group_tags == g,:),...
            nanmean(parameter_matrix(experiment_tags == e & group_tags == 1,:)));
        
        eucdists(1:length(scrap), counter) = scrap;
        
        counter = counter+1;
        
    end
end
%% PCA

[coeff,score,~,~,explained,~] = pca(parameter_matrix); % pca
% Note: By default Matlab Centers data for PCA by subtracting
% the mean from each column (as the means are not quite zero this seems
% appropriate)
[knee_dim] = knee_pt(explained); % Choose this many dimensions

%% Tsne Figure

parameter_matrix_tsne = tsne(parameter_matrix,...
    'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',knee_dim,...
    'Perplexity',30,'Standardize',1,'Verbose',1);

%% FK Export

save_path = '~/.../scn1lab_multidimension/';

csvwrite(strcat (save_path, '/', 'fingerprint.csv'), parametermatrix_fingerprint);
csvwrite(strcat (save_path, '/', 'scatterHom.csv'), parametermatrix_scatterHom);
csvwrite(strcat (save_path, '/', 'eucdists.csv'), eucdists);
csvwrite(strcat (save_path, '/', 'grptags.csv'), group_tags);
csvwrite(strcat (save_path, '/', 'exptags.csv'), experiment_tags);
