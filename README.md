# Simulations for various Matroid Secretary algorithms

The following project consists of simulations for various Secretary and Matroid Secretary algorithms found in literature. 
This is meant to be a framework for testing and exploring further research ideas in this area. Currently, three algorithms have been implemented and are up and working.

1) The classic Secretary problem, described in "Sum the odds to one and stop" by Bruss, F. Thomas [[10.1214/aop/1019160340](https://www.researchgate.net/publication/38351353_Sum_the_odds_to_one_and_stop)]

2) The Secretary problem with stochastic arrivals and departures, as described in the paper "How to Hire Secretaries 
with Stochastic Departures" by Thomas Kesselheim, Alexandros Psomas, Shai Vardi [[arXiV:1909.08660](https://arxiv.org/abs/1909.08660)]

3) An implementation of the Multiple Secretary problem, using the following algorithm: Sample for a fraction of the items arriving, and then choose anything that is in the top - k of the items that have arrived uptil now. This is a natural analogue of "An Optimal Online Algorithm for Weighted Bipartite Matching and Extensions to Combinatorial Auctions" by Thomas Kesselheim, Klaus Radke, Andreas Tönnis and Berthold Vöcking [[10.1007/978-3-642-40450-4_50](http://www.thomas-kesselheim.de/science/onlinematching.pdf)]. The difference is that the sample amount will change based on the value of n and k.

4) The algorithm described in "An Optimal Online Algorithm for Weighted Bipartite Matching and Extensions to Combinatorial Auctions" by Thomas Kesselheim, Klaus Radke, Andreas Tönnis and Berthold Vöcking [[10.1007/978-3-642-40450-4_50](http://www.thomas-kesselheim.de/science/onlinematching.pdf)] has been implemented in the Weighted bipartite matching file. Note that this is actually generalization of the Transversal Matroid Secretary problem, but has the same expected value (1/e, which is optimal).

5) An algorithm for Graphic matroids has been implemented, with good results.

6) An algorithm for Transversal matroids (which crucially depends on the vertex weighted-ness of the problem) has been implemented, with good results.

All of the above simulations should be used to verify that an algorithm is reasonable, and not as a proof that some algorithm is good, as these sample from a distribution and thus are closer to the "Prophet Secretary" problem rather than vanilla Matroid Secretary. To truly verify an algorithm is good using simulations, all kinds of adverserial examples must be tested (which will differ based on the algorithm proposed).
