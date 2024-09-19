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

run_043024:
- added differential fitness weighting for gene segments

run_051324:
- added genotype tracking

run_061324:
- removed tight bottlenecking from within-host replication
- test_sampling.R to find bottleneck size that minimizes drift

run_080224:
- now using gillespie stochastic simulations
- replication is on a per-segment basis instead of a per-genotype basis
- replication propensity takes into account:
  -  likelihood of meeting a matched segment vs a mismatched segment
  -  segment activity
  -  total # of the segment in question



