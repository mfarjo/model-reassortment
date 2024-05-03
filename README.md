# model-reassortment
run_040224:
- gen matrix is randomly assigned from a fitness distribution so original genotypes aren’t balanced
- plotted mean fitness distributions and trajectories at different generations

run_040324:
- starting genotypes are now balanced so that individual gene activities are balanced within a small distribution 
- plotted mean/max fitness distributions and trajectories at different generations

run_041524:
- added mu = 0
- subpopulations can “spread” to take over another subpopulation according to their size 
- added histogram plots for visualization

run_042524:
- took out subpopulation spreading
- initial populations are now in pairs of 2
- each simulation starts with a different pair of genotypes
- added population size metrics
