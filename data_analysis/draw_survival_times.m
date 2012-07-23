function [ ] = draw_survival_times( nb_files, file_name, file_extension )
%DRAW_SURVIVAL_TIMES Draw the survival time repatition of a controler.
% 
%   Default values:
%   * nb_files = 1
%   * file_name = all_survival_tasks_scores
%   * file_extension = dat
%
%   draw_survival_times(nb_files) draws the stats for nb_files# runs
%
%   draw_survival_times(nb_files, file_name) draws the stats for nb_files#,
%                       the data files's names are $file_name$.dat
% 
%   draw_survival_times(nb_files, file_name, file_extension) draws the 
%                       stats for nb_files#, the data files's names are 
%                       $file_name$.$file_extension$
%
%   Examples: 
%      * draw_survival_times()
%      * draw_survival_times(2)
%      * draw_survival_times(1, 'all_survival_tasks_scores_test')
%
%   Dependencies: None
%
%   Expected format of the data file:
%   Generation# Individual# Objective#1 Objective#2 ScoreTask#1 ScoreTask#2... ScScoreTask#n

%% Argument parsing
if nargin < 1
    nb_files = 1;
end

if nargin < 2
    file_name = 'all_survival_tasks_scores';
end

if nargin < 3
    file_extension = 'dat';
end

% If needed, we could add this variables as function parameters
nb_objectives = 2;
nb_lives = 5;


%% Stats computations
survival_tasks_scores = []
for i=1:nb_files
    if (i > 1)
        survival_tasks_scores_temp = load(strcat(file_name, int2str(i), '.', file_extension));
    else
        survival_tasks_scores_temp = load(strcat(file_name, '.', file_extension));
    end
    
    % We only want the first generation
    survival_tasks_scores_temp = survival_tasks_scores_temp(survival_tasks_scores_temp(:, 1) == 0, :);

    survival_tasks_scores_temp = survival_tasks_scores_temp(:, (2 + nb_objectives + 1):(2 + nb_objectives + 1 + nb_lives - 1));
    survival_tasks_scores_temp = reshape(survival_tasks_scores_temp, [], 1);
    
     survival_tasks_scores = horzcat(survival_tasks_scores, survival_tasks_scores_temp);
end


%% Drawing
figure;
hb = boxplot(survival_tasks_scores);
set(gca,'xtick',1:nb_files, 'xticklabel',{'Random','WTA'}) ;
%set(gca,'YTick',[1], {'survival time'}) ;

end

