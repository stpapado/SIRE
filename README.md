# SIRE (Software for IBS and Radiation Effects)
**Description**\
The SIRE code uses the classical Rutherford cross section which is closer to the Piwinski
formalism, for the IBS calculations. It needs as an input the Twiss functions at different locations
of the lattice in order to determine the trajectories of the particles in phase space. This is done
in terms of the two Courant-Snyder and longitudinal invariants, and the 3 phases (betatron and
synchrotron), instead of using the 6 coordinates for position and momentum. The 3 invariants
are conserved between points around the lattice and can only be changed by the effects of IBS,
SR and QE, while the phases are chosen randomly at each given point of the lattice. The time
steps for which the IBS and radiation effects are called should be specified such that they are
larger than the revolution time and smaller than the damping/growth times. Dividing the total
time by the time steps shows how frequently the IBS, SR and QE routines are called. Currently,
the output file giving the evolution of the emittance in all planes presents the values computed
for the specified time steps.

The algorithm SIRE uses to calculate IBS is similar to the one implemented in MOCAC,
where the beam is represented by a large number of macro-particles occupying points in the
6-dimensional phase space. The default distribution defined in SIRE by using a random number
generator, is the Gaussian and is given in action angle variables. So, in order to get a Gaussian
distribution in terms of beam size, the histogram of the macroparticles action angle variables
should form an exponential. In order to apply a different
distribution, someone should either make SIRE generate the proper random deviates or provide
as an input file the action angle variables of all macroparticles for the desired distribution. After
specifying the total beam population and the number of macro-particles, the initial distribution
of the macro-particles can be tracked. The particle distribution in all planes can be saved as
often as requested during the simulation time.\
The steps followed for the IBS effect simulation can be summarized as:\
– For each lattice point defined in the Twiss file, the 3 phases of each macro-particle are
randomly chosen and position and momentum of the macro-particles are calculated.\
– The beam is geometrically divided into a number of cells that is specified for each plane. The
macro-particles are assigned to each cell according to their geometrical position.\
– Based on the classical Rutherford cross section, intra-beam collisions between pairs of
macro-particles are calculated in each cell. The momentum of particles is changed because
of scattering. According to the available computational time, the number of macro-particles
and cells, i.e. the number of collisions each macro-particle experiences, is chosen. The
scattering angles for each collision are determined.\
– The beam distribution is then updated based on the new invariants of the macro-particles.\
– The simulation proceeds to the next lattice point and continues until the end time is reached.

**Adding SR**\
Depending on the elapsed time, the synchrotron radiation damping (RD) acts on the invariants of the macro-particles as an exponential decrement. The routine introduced for this reason
is called after the calculation of the IBS effect for each iteration. Using small iteration time
steps dt (which are much smaller than the damping times and for which the emittances change
adiabatically), the evolution of the transverse emittance and energy spread due to the effects of
IBS and SR can be obtained.

**Lattice recurrences**\
A lattice compression technique named “lattice recurrences”, has been implemented to speed
up the calculations. Since the increase of the invariants due to IBS is linear to the first order
in the traveling time along an element, elements of the full lattice with optics functions differing
less than a specified precision value are considered equal. For such a group of elements, the IBS
effect is evaluated only for one of these elements, resulting in a smaller computational time.

**Convergence studies**\
For a specified set of input beam parameters, various scans should be performed for different combinations of number of macro-particles and cells in order to find the optimal values which provide
a fast tracking and at the same time, guarantee that the scattering process leads to accurate results. In these terms, in order to avoid having a very small number of macro-particles per cell,
the total number of cells is calculated based on the optimal minimum number of macro-particles per cell. A scanning of the total number of cells is performed for an example set of beam parameters
to be used as an input for tracking. The value of the number of macro-
particles per cell after which the variation of the emittances in both planes remains constant is
chosen as the optimal minimum value. After specifying this value, a scanning is performed for a
fixed number of macro-particles, in order to choose the number of cells to be used

**Reduced lattice for LHC**\
As mentioned earlier, one of the inputs required by SIRE are the optical functions along the ring.
As the LHC is a very long accelerator of about 27 km, with a very large number of elements in
the sequence (more than 11000), SIRE requires an extremely long computational time to track
the distribution for all the elements along the ring. Aiming to reduce the computational time,
a study was first performed in order to identify the optimal minimum number of critical IBS
kicks around the lattice, without affecting the overall effect. The IBS growth rates were firstly
calculated for the full optics of the LHC, using the IBS module of the Methodical Accelerator
Design code (MAD-X) which is based on the Bjorken-Mtingwa formalism.
Taking into account the strong IBS kicks along the ring, various lattices with a reduced
number of elements, including the case of the smooth lattice approximation, were tested. Then,
using the IBS module of MAD-X, the emittance evolution was calculated for several sets of beam
parameters to assure that the choice of the elements is valid both for regimes that the effect is
weak and strong. Finally, the optimal lattice chosen consists of only 92 elements for the LHC.
***
**Benchmarking of SIRE with the B-M IBS theoretical model and Data for the LHC**\
The comparison of the code results with analytical formulas has been studied for the
nominal collision energy (7 TeV) for various initial parameters cases. Also, a comparison of
LHC data from Run 2 with simulations performed with SIRE is discussed in. A benchmarking of the Bjorken-Mtingwa (B-M) IBS theoretical model with the SIRE code for
both injection and collision energies is presented for the nominal Batch Compression Merging
and Splitting (BCMS) and the high luminosity (HL-LHC) parameters. For the case
of collision energy, an example showing the comparison between experimental data coming from
Run 2, the SIRE and the B-M analytical formalism, gave excellent results.
***
**Run the code**
See the readme_run for details\
Steps to follow:
1. Compile the code:
  g++ lhcsire.c -o code
2. Run the code :
At least 3 arguments are required. If less then the code gives an error. If 4 arguments then the 4rth one should be the input distribution file.\
  arg1 --> The madx twiss file. Columns has to follow the ordering: name,s1,len1,betx1,alphax1,mux1,bety1,alphay1,muy1,dx1,dpx1,dy1,dpy1\
  arg2 --> The input parameters file. \
  arg3 --> The naming of the output files\
  arg4 --> The input distribution (if the distribution is given as input, if not it assums gaussian), that has a length= # macro-particles, in action angle variables (i.e. in emittance, not in beam size)\
To run:\
  ./code twissfile.tfs paramsfile.dat name\
  or\
  ./code twissfile.tfs paramsfile.dat name distribution.txt
3. Check the outut files when using the arg3 for the naming of the output files we get:\
  RES_TWISS_arg3.txt --> The new twiss file after recurrences\
  RES_EMITTANCE_arg3.txt --> The emittance at each point of the lattice after the IBS kicks (for 1-turn calculations only). Four columns (s, exm, ezm, esm)\
  RES_Growth_Rates_arg3.txt --> File with emittances for each time step. (The Growth rates are the zeroed columns, so they are not saved. ). In this file, L{1}=timesteps (so the NIBSruns), L{2},L{3}=the emittances and energyspread=sqrt(L{4}/2).\
  RES_DISTRIB_arg3.txt --> The distribution file, that has a length= # macro-particles. We can also ask for the output distribution. The particle distribution in all planes can be saved as often as requested during the simulation time.
  
***
**Useful links**\
"Impact of Non-Gaussian Beam Profiles in the Performance of Hadron Colliders" paper by S.Papadopoulou: https://arxiv.org/pdf/1806.07317.pdf \
SIRE code for ABP CWG meeting: https://indico.cern.ch/event/647301/contributions/2630198/attachments/1489047/2313796/ABPCWGpres.pdf  \
SIRE code for ABP CWG webpage: https://twiki.cern.ch/twiki/bin/view/ABPComputing/SIRE  \
"Effect of intrabeam scattering and synchrotron radiation damping when reducing transverse
emittances to augment the LHC luminosity " paper: http://cds.cern.ch/record/1240834/files/sLHC-PROJECT-REPORT-0032.pdf?version=1  \
"Simulation of Intrabeam Scattering" presentation by Vivoli and Martini: https://agenda.linearcollider.org/event/4507/contributions/17682/attachments/14276/23411/CLIC_2010_IBS.pdf  \
"Intrabeam Scattering" paper by Martini and F. Antoniou and Y. Papaphilippou: http://inspirehep.net/record/1507570/files/ICFA69_38-59.pdf \
"OPTICS DESIGN OF INTRABEAM SCATTERING DOMINATED
DAMPING RINGS" thesis of F.Antoniou: https://cds.cern.ch/record/1666863/files/CERN-THESIS-2012-368.pdf \
"Modelling and measurements of bunch profiles at
the LHC" IPAC paper by S.Papadopoulou: https://iopscience.iop.org/article/10.1088/1742-6596/874/1/012008/pdf

