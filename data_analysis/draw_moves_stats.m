function [ ] = draw_moves_stats( nb_actions, file_name )
% DRAW_MOVES_STATS Draws move stats for the survival task.
%
%   Default values:
%   * file_name = 'all_phenotypes_moves.dat' (optional)
%   * nb_actions = 5
%
%   Examples: 
%      * draw_survival_times()
%      * draw_survival_times(2)
%      * draw_survival_times(1, 'all_survival_tasks_scores_test')
%
%   Dependencies: clean_moves.py (MUST BE in the same folder as draw_moves_stats.m
%
%   Expected format of the data file:
%   Generation# Individual# Objective#1 Objective#2 ScoreTask#1 ScoreTask#2... ScScoreTask#n
%   
%   Note: clean_moves.p (called at the beginning of this script) cleans
%   the all_phenotypes_moves.dat log file.


%% Argument parsing
if nargin < 1
    nb_actions = 5;
end

if nargin < 2
    file_name = 'all_phenotypes_moves';
end


% Clean and load the file output by Sferes2
% Make sure that clean_moves.py is in the same folder as draw_moves_stats.m
!python clean_moves.py all_phenotypes_moves
scores = load(strcat(file_name, 'clean.dat'));

%% Stats computations
generations = scores(:, 1);
selected_actions = scores(:, 4);

%Compute the pie
selected_actions_total = [];
for i=1:nb_actions
    selected_actions_total(i) = size(selected_actions(selected_actions(:, 1) == i-1, :), 1);
end

% Compute durations
durations = [];
temp_selected_actions = 4;
cur_idx = 1;
duration = 0;
for i=1:size(selected_actions, 1)
    if (temp_selected_actions == selected_actions(i))
        duration = duration + 1;
    else
        durations = vertcat(durations, [i temp_selected_actions duration]);
        cur_idx = cur_idx + 1;
        duration = 1;
        temp_selected_actions = selected_actions(i);
    end
end


%% Drawing
% figure;
% plot(generations, selected_actions);

figure;
pie3(selected_actions_total);
figure;
pie3(selected_actions_total, {'Wander','Avoid','Reload on Dark', 'Reload on Light', 'Rest'});

figure;
hb = boxplot(durations(:, 3), durations(:, 2));
set(gca,'xtick',1:nb_actions, 'xticklabel',{'Wander','Avoid','Reload on Dark', 'Reload on Light', 'Rest'}) ;

end