
clear all
clearvars
clc

% run all steps in the state estimator sequence

run_powerflow
display(' Pause, when you are ready hit Enter key')
pause

run_scada
display(' Pause, when you are ready hit Enter key')
pause

run_estimator
