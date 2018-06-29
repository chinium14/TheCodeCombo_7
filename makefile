#############################################################################
# 	(C)  	Nitin Rathore and Thomas Knotts: 19 May 2009
#
# Compiling defines options: 
#     -DLJ	        Lennard-Jones interaction potential.
#     -DCOULOMB         Coulombic interaction potential.
#     -DEWALD           Ewald summations for electrostatics.
#     -DNLIST           Neighbor list implementation.
#     -DNEUTRAL         Uses EEF_1 CHARMM parameters with neutral side chains
#     -DWIN             For compiling code in Windows.
#     -DRDIE            Distance dependent dielectric.
#     -DDOS             Density of states calculation.
#     -MMDOS	        Micro-multicanonical Density of states.
#     -DSASA   	        Implicit Solvent Model.
#     -DSASAREX		Replica exchange for co-solvent.
#     -DMEMBRANE	Simulate membrane as osmolyte co-solvent.
#     -DRFC             Reaction Field Correction.
#     -DSHIFTV          Potential shift in Coulombic interaction.
#     -DSHIFTF	        Force shift in Coulombic interactions.
#     -DPRESSURE        Pressure calculations enabled.
#     -DMOVIE           Save pdb files for animation.
#     -DFCC	        Initialize Lennard-Jones particles on FCC lattice.
#     -DMPI             For running code on parallel processors.
#     -DHESSIAN         For computing the hessian numerically.
#     -DEWALD_DEBUG     For print additonal ewald components.
#     -DSTATUS	        For debugging the status of the code.
#     -DNMA	        For running minimization scheme.				
#     -DPR_NPT          For Parnello Rahman NPT simulaiton.
#     -DPR_NPT_DEBUG    Prints additional PR components of conserved quantity.
#     -DSTYPE	        LJ parameters set by site type id instead of site number.
#     -DS14	        Uses special 1-4 LJ parameters and scale Elec by E14FAC.
#     -DCONFIGT         For computing configurational temperature.
#     -DCTDOS           Configurational Temperature Density of states.  		   
#     -DXEDOS	        Density of states in strain space.
#     -DPHI	        For saving the phi files in PHI folder
#     -DCRD             For saving the crd files in CRD folder		
#     -DSTATS	        For saving data for error analysis in DOS sim
#     -DCLIST	        For using a cell list to make neighbor list
#     -DSMD	        For Steered MD simulation with cantilever.			
#     -DNSMD	        For Steered MD simulation without cantilever.
#     -DFLIM	        To limit the forces to some upper limit
#     -DKONS	        To constrain sites to their original positions
#     -DREST	        To restrain sites with harmoic springs
#     -DROTN            For rotating the coordinates at block c	  
#     -DUMBP	        For restraint in z dir only
#     -DIONC            For ion displacement move
#     -DSLIST           For solvent neighbor list	
#     -DBEAD	        Make sigma and epsilon etc. same
#     -DDPCONTROL        Anisotropic scaling of box lengths	
#     -DSPME  	        Smooth particle mesh ewald.  Also need EWALD.(Link w/ fftw3)
#     -DSUG             To use sugar force field parameters of Ha et al.
#     -DTRR 	        To write the gromacs format trajectory in binary	
#     -DFX_EDOS	        To run density of simulaitons with maximum flux	
#     -DGOLIK	        To run Go model of proteins
#     -DDNA_GOLIK       To run Go model of DNA (Must define GOLIK also)
#     -DGOBT            To use bends and torsion in GOLIK(Must define GOLIK also)
#     -DWALL            To add a wall with a GO model
#     -D&WALL           Type of wall where &=N,R,A,C,H
#     -DSPHERE          Type of wall where the shape is sphere
#     -DWALL_DEBUG      To print angle of tethered protein
#     -DSEED            To take the random number seed from system clock
#     -DTDXEDOS	        2 dimensional expanded ensemble density of states
#     -DCAVITY          To add a spherical cavity with a GO model
#     -DACAVITY         To add a spherical cavity with a GO model:felt throughtput inside
#     -DBCAVITY         To add a spherical cavity with a GO model: felt only outside
#     -DBETAPEP         Use this option to have the phi, psi and theta definition go beta peptide
#     -DPREEQUIL        Use with md or replica exhange to ramp temp from 0.001 to sim.T
#     -DDEBHUCK         Use with DNA_GOLIK for a debye-huckel screening of coulombics
#     -DFSH             Use force shift for coulombics in DNA_GOLIK
#     -DTWHAM           Use to accumulate histograms for temperature replica exchange.
#     -DXWHAM           Uset to accumulate histograms for umbrella sampling
#     -DZHOU		Use to compute per-atom stress tensors
#     -DZEROAM		Use to zero out the angular momentum along with the linear momentum
#     -DCHECKP          Generate and read checkpoint files automatically (useful for submitting jobs with preemptive scheduling)
#     -DOHEADER		Output .header files to indicate what the fields in the other output files mean
#     -DSCRIPT		For reading the new format of simul.input.
#==========================================================================================

#-------------------------------------------------------------------
# Select a compiler and options to use (comment appropriate lines)
#-------------------------------------------------------------------
#For parallel compiling
#<<<<<<< makefile
CC       =  gcc
#CC	=  mpicc
#CC       =  pgcc #don't forget to put -DMPI in DFLAGS
#CFLAGS   = -fast -ta=nvidia,cc20 -Minfo -c 
CFLAGS   = -O3 -ffast-math -funroll-all-loops -Wall -std=c99  
#CFLAGS   = -O2 -ffast-math -finline-functions -funroll-all-loops -Wall -std=c99 -g -pg

#=======
#CC       =  mpicc #don't forget to put -DMPI in DFLAGS
#CFLAGS   = -fast -w -static-libcxa -xW -axW -std=c99 
#CFLAGS   = -finline-functions -funroll-all-loops -std=c99 -static

#CC       =  gcc #don't forget to put -DMPI in DFLAGS
#CFLAGS   = -O3 -ffast-math -finline-functions -funroll-all-loops -Wall -std=c99
#>>>>>>> 1.5
#Using gcc
#<<<<<<< makefile
#CC       =  gcc
#CFLAGS   = -O3 -ffast-math -finline-functions -funroll-all-loops -Wall -std=c99 
#=======
#CC       =  gcc
#CFLAGS   = -g -std=c99 -O3 -ffast-math -finline-functions -funroll-all-loops -Wall -std=c99 

#Debugging with gcc
#CC       =  gcc
#CFLAGS   =  -Wall -std=c99 -g
#>>>>>>> 1.5

#Using icc 
#CC	  = icc #Can also define -DMKL and CHANGE LIBRARY -lm TO -lmkl_core
#CFLAGS   = -fast -w -static-libcxa -xW -axW -std=c99
#CFLAGS    = -w -O2 -static-libcxa -std=c99 

#------------------------------------------------------------------
# Select the flags to define
#------------------------------------------------------------------
#Umbrella sampling, Brooks Go model with a sphere
#DFLAGS    = -DNLIST -DDLIST\
            -DBGO -DSEED \
            -DREST -DWALL -DSPHERE\
            -DTRR -DKONS -DXWHAM

#Sphere 
#DFLAGS    = -DNLIST -DDLIST\
            -DBGO -DSEED -DTWHAM\
            -DTRR -DWALL\
            -DSPHERE  # -DZEROAM #-DMPI
#Sphere replica tethered
#DFLAGS    = -DNLIST -DDLIST\
            -DBGO -DSEED -DTWHAM\
            -DTRR -DREST -DWALL\
            -DSPHERE  -DMPI # -DZEROAM #-DMPI

# Umbrella sampling, Brooks Go model with a wwall
#DFLAGS    = -DNLIST -DDLIST \
            -DBGO -DSEED -DXWHAM -DTRR\
            -DREST -DWALL -DWWALL -DKONS 
#DFLAGS    = -DNLIST -DDLIST \
            -DBGO -DSEED -DXWHAM -DTRR\
            -DREST 

# Replica exchange, Brooks Go model without a wwall

# RE with different solvent conditions.
DFLAGS    = -DNLIST -DDLIST\
            -DBGO -DSEED -DTRR\
            -DWALL -DWWALL
#-DTWHAM -DMPI

#DFLAGS    = -DNLIST -DDLIST \
            -DBGO -DSEED -DTRR\
            -DSASA -DMEMBRANE \
            -DTWHAM -DMPI#-DMEMINSIDE #-DMPI

#DFLAGS    = -DNLIST -DDLIST \
            -DBGO -DSEED -DTRR \
            -DSASA -DMEMBRANE #-DMEMINSIDE

#DFLAGS    = -DNLIST -DDLIST \
            -DBGO -DSEED -DTRR \
            -DREST -DSASA -DWALL -DWWALL

#DFLAGS    = -DNLIST -DDLIST\
            -DBGO -DSEED -DTRR\
            -DSASA -DMPI -DTWHAM

# MD, CHARMM19
#<<<<<<< makefile
#DFLAGS    = -DNLIST -DDLIST -DTLATE -DTRR -DCOULOMB -DSASA -DNEUTRAL -DRDIE\
            -DSHIFTF -DMKL -DMPI
#=======
#DFLAGS    = -DNLIST -DDLIST -DNEUTRAL -DCOULOMB -DSASA -DRDIE -DSHIFTF -DTLATE -DTRR 
            
#>>>>>>> 1.5

#MD, DNA model of Knotts et. al. using Debye Huckel for coulombics


#DFLAGS    = -DNLIST -DDLIST -DSTYPE -DS14 -DTRR -DCOULOMB -DRDIE -DNEUTRAL\
	    -DSCRIPT -DZHOU -DOHEADER -DCHECKP -DPRESSURE -DSASA -DSHIFTF


#DFLAGS    = -DNLIST -DDLIST -DMPI -DTWHAM -DBGO -DPRESSURE -DSCRIPT -DCHECKP -DZHOU -DOHEADER

#MD, CHARMM22,27, Smooth Particle Mesh Ewald
#DFLAGS   =	-DNLIST -DDLIST -DTWHAM -DBGO -DMPI -DSEED -DTRR 
#DFLAGS   =	-DNLIST -DDLIST -DBGO -DSEED -DTRR 

#<<<<<<< makefile
#DFLAGS = -DNLIST -DDLIST -DTRR -DSEED -DPRESSURE -DLJ  
#=======
#DFLAGS = -DNLIST -DDLIST -DTRR -DSEED -DPRESSURE -DLJ -DFCC 
#>>>>>>> 1.5
 
# --------------------------------------------------------------------------------
# Select include and library files
# --------------------------------------------------------------------------------

# For Debye
#<<<<<<< makefile
#INCL      = -I/usr/include/ -I/share/apps/include/ -I/share/apps/mpich2-1.0.8/lib \
            -I/share/apps/intel/mkl/10.0.4.023/include
#=======
#INCL      = -I/usr/include/ -I/share/apps/include/  -I/share/apps/mpich2-1.0.8/lib \
            -I/share/apps/intel/mkl/10.0.4.023/include
#>>>>>>> 1.5

#<<<<<<< makefile
#LIBS    = -L/usr/lib64  -L/share/apps/lib/ -L/share/apps/mpich2-1.0.8/lib \
          -L/share/apps/intel/mkl/10.0.4.023/lib/em64t -lmkl_em64t -lfftw3 #-lefence
#=======
#LIBS    = -L/usr/lib64  -L/share/apps/lib/ \
	   -L/share/apps/intel/mkl/10.0.4.023/lib/em64t -lm -lfftw3 # -lefence
#>>>>>>> 1.5

# For Marylou
INCL      = -I/usr/include/

LIBS    = -L/usr/lib64 -lm #-lfftw3 

# -----------------------------------------------------------------------------
# Name of the executable
# -----------------------------------------------------------------------------

#<<<<<<< makefile
#EXEC	  = md_0.81_0.768_0.20_0.01_0.05
#=======
#EXEC	  = replica_sphere_400A
#EXEC	  = umb_sphere_400A_new
#EXEC	  = sphere_80n
#EXEC 	  = replica_sphere_20A_scaled
#EXEC 	  = replica_sphere_400A_new
#EXEC 	  = newwall_hybrid
EXEC	 = wall_new
#EXEC	 = replica_solvent1011
#EXEC	 = replica_solvent111_wall
#EXEC	 = md_solvent_wall451_test
#EXEC	 = md
#EXEC	 = md_membrane_inside
#EXEC	 = md_membrane
#EXEC	 = replica_membrane
#EXEC	 = replica_membrane_inside


#>>>>>>> 1.5

#
#############################################################################
# nothing should be changed below here


SRCS =  allocate.c avrg_xe_mpi.c                                                 \
        brent.c boxinv.c boxscale.c                                              \
        calcvalue.c cbend.c cbond.c cnonbond_nblist.c ctorsion.c csasa.c         \
         ctdos.c ctdos_hmc.c ctdos_pivot.c ctdos_trans.c cbend2.c cbond2.c       \
         cimproper2.c csasa2.c csasab_bgo.c cnonbond_nblist2.c ctorsion2.c configtemp.c       \
         checks.c crestraint.c clist.c csmd.c curey.c curey2.c                   \
         cimproper.c cplx.c curr_status.c cgolik.c cwall.c ccavity.c             \
         crestraint2.c cgolik2.c cwall2.c ctdos_assoc_trans.c checkpoint.c       \
        dgauss.c drift.c dos_awrite.c dos_wwrite.c dos_vblock.c                  \
         dos_svalues.c dos_stats.c dos.c dos_interpolate.c dos_assoc_trans.c     \
        ebond_mc.c ebend_mc.c etors_mc.c eurey_mc.c eimpr_mc.c erest_mc.c        \
         esasa_mc.c enbnd_mc.c ewald.c end_to_end.c ewald2.c ewald_mc.c          \
         ebond2_mc.c ebend2_mc.c etors2_mc.c enbnd2_mc.c eimpr2_mc.c             \
         esasa2_mc.c ewald2_mc.c egolik_mc.c egolik2_mc.c erest2_mc.c            \
         eurey2_mc.c                                                             \
        forces.c force_long.c force_short.c f1dim.c frprmn.c forces_mc.c         \
         fx_edos.c fx_edos_hmc.c fx_edos_pivot.c fx_edos_trans.c                 \
        hybrid_MC.c helix.c hessian.c                                            \
        integrate_vverlet.c init.c integrate.c integrate_nhc.c isokin.c          \
         integrate_nhcp_full.c integrate_npt_full.c integrate_smd.c              \
         integrate_dos_mts.c integrate_mts.c integrate_mts_nhc.c                 \
         integrate_mts_npt_full.c ioflush.c                                      \
        kinet.c                                                                  \
        linmin.c ljpot.c                                                         \
        main.c mmdos_dke.c mmdos.c mmdos_hmc.c minim.c movie.c                   \
         mmdos_pivot.c min_nma.c mnbrak.c matrixmath.c min_image.c               \
         mcmove_bkp.c mcmove_trans.c mcmove_pivot.c mcmove_displ.c               \
         mcmove_hmc.c mem_free.c                                                 \
        nblist.c npt_md.c nve_md.c nvt_md.c nblist_pivot.c native_cont.c         \
         nvt_mc.c                                                                \
        ofile.c outend.c output.c output_headers.c                               \
        pbc_all.c pbc.c pbc_init.c pivot.c pbc_chk.c pme_setup.c                 \
         pme_calc.c pme_qgrid.c pme_spline.c preequil.c                          \
        quatfit.c                                                                \
        random.c read_atom.c read_config.c read_topology.c read_sim.c            \
         read_veloc.c repexch.c rad_gyr.c read_cons.c read_rest.c                \
         rotate_xyz.c ran_int.c read_golik.c read_wall.c read_swap.c             \
         read_dnagolik.c                                                         \
        save_config.c save_veloc.c setlj.c svalues.c swap_dos.c                  \
         swap_mm_mpi.c swap_ct_dos.c swap_ct_mpi.c swap_xe_dos.c                 \
         swap_xe_mpi.c smd_work.c swap_box.c swap_mpi.c swap_mm_dos.c            \
         swap_box_mpi.c stress.c                                                 \
        tdxedos_axial.c translate.c trans.c tdxedos.c tdxedos_init.c             \
         tdxedos_hmc.c tdxedos_linear.c tdxedos_scale.c tdxedos_solv.c           \
         tdxedos_trans.c tdxedos_pivot.c twham.c tak_histogram.c                 \
        update_weights.c                                                         \
        vblock.c velscale.c vinit.c virialcor.c vinit_hyb.c                      \
        weight.c write_gro.c write_param.c write_trr.c wall_angle.c              \
        xyz_config.c xedos.c xedos_hmc.c xedos_axial.c xedos_trans.c             \
         xedos_scale.c xedos_linear.c xedos_pivot.c xedos_solv.c                 \
         xedos_cosol.c xwham.c pbc_npt_full.c hydro_index_bgo.c			 \
	wall_nu_nu.c wall_sphere.c csasa_bgo.c vdw_radii_bgo.c read_solvent.c
	
     
OBJS = ${SRCS:.c=.o}

all: $(EXEC)

%.o: %.c
	${CC} ${CFLAGS} ${DFLAGS} ${INCL} -c  $< -o $@

$(EXEC):  ${OBJS}
	$(CC) ${CFLAGS} ${DFLAGS} ${INCL} -o $@ ${OBJS} $(LIBS)
	echo $(EXEC)

#$(EXEC):  ${OBJS}
#	$(CC) ${CFLAGS} ${DFLAGS} ${INCL} ${OBJS} $(LIBS)
#	echo $(EXEC)

clean:
	rm -f *.o
	rm -f $(EXEC)

