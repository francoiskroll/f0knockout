% Vp_Extract 

% marcus.ghosh.11@ucl.ac.uk 

% FK FK FK
% just the ViewPoint errors fix module
% i.e. fixes the ViewPoint tracking glitches and export the data

% COMMENTS
% can smooth, but might want to do trial/error (as short timeframes, cannot
% smooth the same way as the long experiments)
% ! change fps below if not 25 fps

% v1: also export delta_px_sq, better if want to bin into a barplot

%% Info 

% Extracts data from Viewpoint files with data for every frame

% Input
    % Excel Sheets - output from re-running a viewpoint file  
        % Each with:
        % Rows - every frame for every animal
        % Columns including:
            % Type - data (101) or errors (~=101)
            % Location - which animal each data point came from
            % Data - pixels changed that frame
     % Genotype list (txt): 
        % Columns - group names  
        % Rows - ROI ID's assigning each fish to a group 
    
%% Assumptions 

% Data comes from 192 regions of interest (ROIS)

%% Dependencies 

% Files  
    % dir2
        % https://github.com/ghoshm/Utilities/blob/master/dir2.m 

    % lbmap
        %https://uk.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps

    % Nat Sort Files 
        %http://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort 

    % ProgressBar
        %http://uk.mathworks.com/matlabcentral/fileexchange/6922-progressbar
    
% MATLAB Products 
    % Statistics and Machine Learning Toolbox
    % Curve Fitting Toolbox
    % Parallel Computing Toolbox
    % MATLAB Distributed Computing Server # I installed the Amazon Services
    % one

%% Options
    % For *'s see "Notes on Options" below 
    
% General  
lines_per_sheet = 50000; % Specify the number of data points per Excel sheet 
box = 2; % set which of the two boxes you want to use (*) 
threshold = 200; % Maximum delta px value (**)  
top_up = []; % alter to light boundaries where you topped up fish water.
    % E.g. Day 2 and 3 (top_up = [2 4]). 
    % E.g. Not topped up (top_up = []). 
top_up_threshold = 8; % threshold for cropping data where water was topped up (***)
time_bins = 60; % Choose plot smoothing (from seconds) 
days = [1]; % number of days % typically 1 2 3 4
nights = []; % number of nights % typically 1 2 3

fps = 25

% Colors
col = 'RedBlue'; % (****) 
night_color = [0.9608 0.9608 0.9608]; % color of night shading  

%% Notes on Options  
%* at this point its easier to run each seperatly then merge
    % the data later using Vp_Analyse  
    
% ** A hard threshold is a simple solution, though removing whole bouts
    % including these values is best as it avoids leaving "cut bouts" 
    % From my 3 WT experiments (with lids so no hands), the maximum value I
    % Observe is 165px.
    % From a PTZ Dose response (you may expect higher values) the data
    % doesn't exceed this
    % From looking @ 6 Hcrt experiments and the PTZ experiment - data(data > 0)
    % prctile and std, there is no obvious flexible cut off that could be used. 
    
% *** Again a hard coded threshold is a simple solution. 
    % So far experiments @ 25Hz work with 10 
    % experiments @ 15Hz work with 8

% **** Choice of 
    % 'Blue'       Single-hue progression to purlish-blue (default)
    % 'BlueGray'   Diverging progression from blue to gray
    % 'BrownBlue'  Orange-white-purple diverging scheme
    % 'RedBlue'    Modified spectral scheme
    
%% Selecting Files  

% Select a folder of Excel Sheet
disp ('>Folder with Excel files?')
folder_path = uigetdir([],'Select a folder of Excel sheets'); % Choose your experiment
folder_open = dir2(folder_path); % Open this folder
disp(horzcat('Running File ',folder_path)); % Report file choice  

% Select a geno_list
disp('>Metadata file?')
[filename, pathname] = uigetfile('*.txt', 'Select a Genotype List'); % Select a geno file
if isequal(filename,0) % If no file is selected
    error('No File Selected') % Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) % Show selected filename
end
geno_list = importdata(strcat(pathname,filename),'\t',2); % Load genolist

% Select a save path 
disp ('>Save output?')
save_pathname = uigetdir([],'Select a save location'); 
[~,save_name,~] = fileparts(filename); % assign file name for save 
disp(['Save path ',save_pathname]); % report 

%% Load Data from Excel Sheets 

tic
% Pre-allocation
time = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % time {1}
%pause(30); % Wait for memory 
data_type = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Data type {2} 
%pause(30); % Wait for memory 
fish_id = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Fish id {3}
%pause(30); % Wait for memory 
delta_px = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Delta px {4}
%pause(30); % Wait for memory 
sheet_names = cell(size(folder_open,1),1); % Excel sheet names  

% Ordering by File Name 
for f = 1:size(folder_open,1) % For each excel file 
    sheet_names{f} = folder_open(f).name; % Take it's name 
end % Note that these will be in "computer" order 
    % Ie. 1-10, 100, 1000 etc 

[~,O] = natsortfiles(sheet_names); % Use natsortfiles to sort by file name
    clear sheet_names; % Clear Sheet names 
    
a = 1; % Start a counter 
progress = 0; % Start a timer
data_type_errors = 0; % Start a counter  
order_errors = 0; % Start a counter 
order_errors_size = []; % pre-allocate an empty vector
progressbar('Files') %Initialise progress bars 
for f = O' % For each Excel file
    
    % Load data 
        % Raw data structure 
            % 1 - Time 
            % 2 - Data type 
            % 3 - Fish ID 
            % 4 - Delta px
    fid = fopen(strcat(folder_path,'/',folder_open(f).name)); % Open it % FK ! I changed \ to /
    if f == O(1) % For the first file Skip the first blank lines
        raw_text = textscan(fid,... % Read in the data
            '%*f32 %f32 %*f32 %f32 %*1s %s %f32','headerlines',192 + 1);
    else % For the rest of the files 
        raw_text = textscan(fid,... % Read in the data
            '%*f32 %f32 %*f32 %f32 %*1s %s %f32','headerlines',1);
    end
    
    % Error Handling 
    
    % Dealing with data_type errors (e)
        % Note that this includes cropping the end of the experiment off 
    found = find(raw_text{2} ~= 101); % Find errors
    if isempty(found) ~= 1 % If there are errors 
        for e = 1:size(raw_text,2) % For each data type being imported
            raw_text{e}(found) = []; % Remove errors
        end
        data_type_errors = data_type_errors + 1; % Add to error counter  
    end
    clear found;
    
    % Dealing with ordering errors 
    fish_id_order = str2num(char(raw_text{3})); % Convert str to num
    fish_id_order_check = diff(fish_id_order); % Diff these values
    found = find(fish_id_order_check ~= 1 & fish_id_order_check ~= -191 &...
        isnan(fish_id_order_check) == 0,1,'first'); % Find the first frame drops 
    
    while isempty(found) == 0 % While there are dropped frames 

        % Use a sliding window to find where the pattern normalises
        % Cut backwards
        found_clip_b = [191 1];
        if found - found_clip_b(1) < 1
            found_clip_b = [found + 1 found + 1]; 
        else
            while sum(fish_id_order_check(found - found_clip_b(1):found - found_clip_b(2))) ~= 191
                found_clip_b = found_clip_b + 1; % Slide window
                
                % Catch Running past the start of the file exception
                if found - found_clip_b(1) < 1
                    found_clip_b = [found_clip_b(1) + 1 found_clip_b(1) + 1];
                    break
                end
            end
        end
        
        % Cut forwards
        found_clip_f = [1 191];
        if found + found_clip_f(2) > length(fish_id_order_check)
            found_clip_f = [length(fish_id_order_check) - found + 2 ... 
                length(fish_id_order_check) - found + 2]; 
        else
            while sum(fish_id_order_check(found + found_clip_f(1):found + found_clip_f(2))) ~= 191
                found_clip_f = found_clip_f + 1; % Slide window
                
                % Catch Running past the end of the file exception
                if found + found_clip_f(2) > length(fish_id_order_check)
                    found_clip_f = [found_clip_f(2)+1 found_clip_f(2)+1];
                    break
                end
            end
        end
        
        % Now set values between these sections to NaN
        fish_id_order(found - found_clip_b(2)+2:found + found_clip_f(1)-1) = NaN;
        order_errors_size(order_errors+1,1) = ...
            size(found - found_clip_b(2)+2:found + found_clip_f(1)-1,2); 
            % Store the size of the removed data 
        clear found_clip_b found_clip_f
        
        fish_id_order_check = diff(fish_id_order); % Diff the new fish order
        found = find(fish_id_order_check ~= 1 & fish_id_order_check ~= -191 &...
            isnan(fish_id_order_check) == 0,1,'first'); % Find other frame drops
        order_errors = order_errors + 1; % Add to the order errors counter 
    end
        
    % Data Storage 
    if a == 1 % For the first file 
        time(a:size(raw_text{1},1),1) = raw_text{1}; % Time 
        data_type(a:size(raw_text{1},1),1) = raw_text{2}; % Data Type 
        fish_id(a:size(raw_text{1},1),1) = fish_id_order; % Fish Id order 
        delta_px(a:size(raw_text{1},1),1) = raw_text{4}; % Delta px 
            a = a + size(raw_text{1},1); % Add to counter 
    else % For all other files 
        time(a:a+size(raw_text{1},1)-1,1) = raw_text{1}; % Time 
        data_type(a:a+size(raw_text{1},1)-1,1) = raw_text{2}; % Data Type 
        fish_id(a:a+size(raw_text{1},1)-1,1) = fish_id_order; % Fish Id order 
        delta_px(a:a+size(raw_text{1},1)-1,1) = raw_text{4}; % Delta px 
            a = a + size(raw_text{1},1); % Add to counter
    end
        
    % Clean up 
    fclose(fid); clear fish_id_order fish_id_order_check found raw_text; % Clear variables
    
    progress = progress + 1; % Add to timer 
    progressbar(progress/size(folder_open,1)); %Update the Fish progressbar
    
end
disp('Time taken to load data from Excel sheets'); % report progress
toc 

% Report back on errors 
disp(horzcat(num2str(data_type_errors),' Files with Data Type Errors')); 
disp(horzcat(num2str(order_errors),' Order Errors')); 
for e = 1:size(order_errors_size,1) % For each order error 
    disp(horzcat('Order Error ',num2str(e),' Size = ',num2str(order_errors_size(e)))); 
end 
disp(horzcat('Ran File ',folder_path)); % Report file choice  

clear a e ans data_type data_type_errors f fid folder_open folder_path ...
    lines_per_sheet O order_errors order_errors_size progress

%% Checking for Persistant Ordering Errors
% At this point only one type of ordering error will remain
% This will be the rare case where the ordering error bridges two Excel sheets
% E.g. 171030_16

% Check for remaining Errors
scrap = diff(isnan(fish_id));
% This will leave either
% +1 (A number to a NaN)
% -1 (A NaN to a number)
% 0 (number to number)

% Check for Backwards errors
check(:,1) = find(scrap == 1); % Locs
check(:,2) = fish_id(scrap == 1); % Values (should all be 192)

if isempty(find(check(:,2) ~= 192)) == 0 % If there are errors
    err = find(check(:,2) ~= 192); % define them
    fish_id_diff = diff(fish_id); % diff fish_id
    
    for e = 1:size(err,1) % for each error
        
        % Use a sliding window to find where the pattern normalises
        % Cut backwards
        found = check(err(e,1),1); % set the error location
        found_clip_b = [191 1];
        if found - found_clip_b(1) < 1
            found_clip_b = [found + 1 found + 1];
        else
            while sum(fish_id_diff(found - found_clip_b(1):found - found_clip_b(2))) ~= 191
                found_clip_b = found_clip_b + 1; % Slide window
                
                % Catch Running past the start of the file exception
                if found - found_clip_b(1) < 1
                    found_clip_b = [found_clip_b(1) + 1 found_clip_b(1) + 1];
                    break
                end
            end
        end
        
        fish_id(found - found_clip_b(2)+2:found) = NaN;
    end
    
end

clear check err e found found_clip_b fish_id_diff

% Check for forwards errors
check(:,1) = find(scrap == -1) + 1;
check(:,2) = fish_id(find(scrap == -1)+1);

if isempty(find(check(:,2) ~= 1)) == 0 % If there are errors
    err = find(check(:,2) ~= 1); % define them
    fish_id_diff = diff(fish_id); % diff fish_id
    
    for e = 1:size(err,1) % for each error
        
        % Cut forwards
        found = check(err(e,1),1); % set the error location
        found_clip_f = [1 191];
        if found + found_clip_f(2) > length(fish_id_diff)
            found_clip_f = [length(fish_id_diff) - found + 2 ...
                length(fish_id_diff) - found + 2];
        else
            while sum(fish_id_diff(found + found_clip_f(1):found + found_clip_f(2))) ~= 191
                found_clip_f = found_clip_f + 1; % Slide window
                
                % Catch Running past the end of the file exception
                if found + found_clip_f(2) > length(fish_id_diff)
                    found_clip_f = [found_clip_f(2)+1 found_clip_f(2)+1];
                    break
                end
            end
        end
        
        fish_id(found:found + found_clip_f(1)-1) = NaN;
    end
    
end

clear check err fish_id_diff e found found_clip_f fish_id_diff scrap

disp('Found & Corrected Persistant Ordering Errors'); % Report

%% Reshape The Data 

% Check that the tracking start's with fish 1 
if fish_id(1) ~= 1 % if not 
   crop = find(fish_id == 1,1,'first') - 1; % find just before fish 1   
   % Crop all of the data 
       delta_px(1:crop) = []; 
       fish_id(1:crop) = []; 
       time(1:crop) = []; 
   clear crop; 
end 

% Check the number of frames per fish is correct 
frames_per_fish = zeros(1,max(fish_id)); 
for f = 1:max(fish_id) % For each fish 
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    frames_per_fish(f) = size(found,1); % Store number of frames 
    disp(horzcat('Found frames for fish Number ',num2str(f),' of ',...
        num2str(max(fish_id)))); % Report on progress 
end 

if min(frames_per_fish) ~= max(frames_per_fish) % If the number of frames is not equal 
   error('Data formatted Incorrectly'); % Call an error  
end 

% Delta px sq 
for f = 1:max(fish_id) % For each fish  
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    
    if f == 1 % For the first fish 
        delta_px_sq = nan(frames_per_fish(1),max(fish_id),'single'); 
        % Pre-allocate 
    end
    
    delta_px_sq(:,f) = delta_px(found); % Take delta_px values
    disp(horzcat('Collected frames for fish ',num2str(f),...
        ' of ',num2str(max(fish_id)))); % Report on progress 
end 
clear f found delta_px

% Time sq 
    % Note this is separated from delta px sq for ease of memory handling 
for f = 1:max(fish_id) % For each fish  
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    
    if f == 1 % For the first fish 
        time_sq = nan(frames_per_fish(1),max(fish_id),'single'); 
        % Pre-allocate 
    end
    
    time_sq(:,f) = time(found); % Take time values 
    disp(horzcat('Collected time for fish ',num2str(f),...
        ' of ',num2str(max(fish_id)))); % Report on progress 
end
clear f found time fish_id frames_per_fish

%% Organise the data 

if box == 1
    delta_px_sq(:,97:end) = []; % Remove the unused box
else
    delta_px_sq(:,1:96) = []; % Remove the unused box  
end

delta_px_sq = delta_px_sq - 1; % Set minimum value to zero 

%% Remove "Noise" - set values to zero
% 1. Abnormally high viewpoint values
% 2. Topping up Fish Water

% 1. Remove High Viewpoint values
% Note that this is adapted from the parameter extraction code (below)
wake_cells = cell(1,size(delta_px_sq,2)); % Wake Cells (bout parameters)

% Finding transitions
delta_px_sq_scrap = delta_px_sq;
delta_px_sq_scrap(delta_px_sq_scrap > 0) = 1; % Find active frames
delta_px_sq_scrap = diff(delta_px_sq_scrap); % Diff to find transitions
% 1 = inactive to active
% -1 = active to inactive

for f = 1:size(delta_px_sq,2) % For each fish
    % Note this this runs apporximately twice as fast as using a
    % For loop
    
    % Starts - ensures no bouts are lost at the start
    if  delta_px_sq(1,f) > 0 % If active in first bin
        wake_cells{1,f}(:,1) = [1 ; find(delta_px_sq_scrap(:,f) == 1)+1]; % Find active bout starts
    else % Ie. if inactive in first bin
        wake_cells{1,f}(:,1) = find(delta_px_sq_scrap(:,f) == 1)+1; % Find active bout starts
    end
    
    % Ends - ensures no bouts are lost at the end
    if delta_px_sq(size(delta_px_sq,1),f) > 0 % If active in last bin
        wake_cells{1,f}(:,2) = [find(delta_px_sq_scrap(:,f) == - 1);...
            size(delta_px_sq,1)]; % Find active bout ends
    else
        wake_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == - 1);
    end
    
    % Parameter extraction
    wake_cells{1,f}(:,3) = NaN; % Pre-allocate
    
    % Active bouts
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        wake_cells{1,f}(b,3) = nanmax(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Max
        
        if wake_cells{1,f}(b,3) > threshold % Hard coded threshold
            delta_px_sq(wake_cells{1,f}(b,1):...
                wake_cells{1,f}(b,2),f) = 0; % Bin to zero
        end
        
    end
    
end

clear b delta_px_sq_scrap f wake_cells threshold

% 2. Filter out Hands & Truncated Bouts - V2 
%figure; plot(nanmax(delta_px_sq')); title('Max'); 
%figure; plot(nanmean(delta_px_sq')); title('Mean'); 

if isempty(top_up) == 0 % if fish h20 was topped up
    top_up_bin = nan(size(top_up,2),2,'single'); % pre-allocate (top ups x start/stop)
    
    for t = 1:size(top_up,2) % For each top up
        [~,top_up_bin(t,1)] = find(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)'))) >= ...
            nanmean(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))) + ...
            top_up_threshold*nanstd(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))),1,'first'); % find start 
        
        [~,top_up_bin(t,2)] = find(abs(diff(max(delta_px_sq(lb(top_up(t)):(lb(top_up(t)) + top_up_bin(t,1) + (fps*60*20)),:)'))) >= ...
            nanmean(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))) + ...
            top_up_threshold*nanstd(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)')))),1,'last'); % find stop 
        % Note that here I only look within 20mins (this helps avoid
        % reminaing viewpoint glitches "(fps*60*20)") 
        
        figure; hold on;
        plot(abs(diff(max(delta_px_sq(lb(top_up(t)):lb(top_up(t)+1),:)'))));
        plot([top_up_bin(t,1) top_up_bin(t,1)] - (fps*90),[0 200],'r','linewidth',3);
        plot([top_up_bin(t,2) top_up_bin(t,2)] + (fps*90),[0 200],'r','linewidth',3);
        
        % To ensure you get all of the noise, cut a bit more either side 
        top_up_bin(t,1) = lb(top_up(t)) + top_up_bin(t,1) - (fps*90); % go 90s further back
        top_up_bin(t,2) = lb(top_up(t)) + top_up_bin(t,2) + (fps*90); % go 90s further forwards
        
        for f = 1:size(delta_px_sq,2) % for each fish
            if delta_px_sq(top_up_bin(t,1),f) == 0 && delta_px_sq(top_up_bin(t,2),f) == 0
                delta_px_sq(top_up_bin(t,1):top_up_bin(t,2),f) = 0; % set these values to zero
            else % if they have bouts overlapping with these cuts
                delta_px_sq(top_up_bin(t,1)- ...
                    (find(flip(delta_px_sq(1:top_up_bin(t,1),f)) == 0,1,'first')-2):...
                    top_up_bin(t,2) + (find(delta_px_sq(top_up_bin(t,2):end,f) == 0,1,'first')-2),f) = 0;
            end
        end
        

        
    end
    
else % if not topped up
    top_up_bin = []; % store a blank variable
end

clear f t top_up

%figure; plot(nanmax(delta_px_sq')); title('Max - Post '); 
%figure; plot(nanmean(delta_px_sq')); title('Mean - Post '); 

%% Group the data by condition 

% Generate group tags 
group_tags = nan(size(delta_px_sq,2),1); % Pre-allocate
for g = 1:size(geno_list.data,2) % For each group 
    group_tags(geno_list.data(1:find(isnan(geno_list.data(:,g))==0,1,'last'),...
        g)) = g; % Assign group membership  
end 
delta_px_sq(:,isnan(group_tags)) = []; % Remove data   
group_tags(isnan(group_tags)) = []; % Remove blank values 

clear g  

%% Smoothing data into seconds
% Pre-allocate
delta_px_sq_sec = nan(size(1:(fps):size(delta_px_sq,1),2),...
    size(delta_px_sq,2),'single'); % time (seconds) x fish
delta_px_sq_sec_smooth = nan(size(1:(fps):size(delta_px_sq,1),2),...
    size(delta_px_sq,2),'single'); % time (seconds) x fish

% Smooth each fish's activity & parameters
for f = 1:size(delta_px_sq,2) % For each fish
    
    a = 1; % Start a counter
    for t = 1:(fps):size(delta_px_sq,1) % For each second
        if t + (fps-1) < size(delta_px_sq,1) % Check to prevent running off the end
            delta_px_sq_sec(a,f) = nansum(delta_px_sq(t:t+(fps-1),f)); % Bin activity
        
        a = a + 1; % Add to counter
        
        else 
            delta_px_sq_sec(a,f) = 0; % This prevents smoothing errors
        end
        
    end
    
    % Smooth the activity data
    delta_px_sq_sec_smooth(:,f) = smooth(delta_px_sq_sec(:,f),time_bins);
    disp(horzcat('Smoothed fish ',num2str(f),' of ' ,...
        num2str(size(delta_px_sq,2)))); % Report progress
end

clear a f p t 

%% Export delta_px

save_name2 = strrep(save_name, 'genotype', '');

csvwrite (strcat (save_pathname, '/', save_name2, '_deltapxsq.csv'), delta_px_sq);
csvwrite (strcat (save_pathname, '/', save_name2, '_deltapxsqsec.csv'), delta_px_sq_sec);
csvwrite (strcat (save_pathname, '/', save_name2, '_deltapxsqsecsmooth.csv'), delta_px_sq_sec_smooth);

%% Output from Vp_Extract to csv
% FK 12/03/2019 

% save_name2 = strrep(save_name, '_DATA', '');

% csvwrite (strcat (save_pathname, '/', save_name2, '_deltapxsecsmooth.csv'), delta_px_sq_sec_smooth)
% csvwrite (strcat (save_pathname, '/', save_name2, '_lbsec.csv'), lb_sec)

% csvwrite (strcat (save_pathname, '/', save_name2, '_deltapx.csv'), delta_px_sq)
% csvwrite (strcat (save_pathname, '/', save_name2, '_lb.csv'), lb)