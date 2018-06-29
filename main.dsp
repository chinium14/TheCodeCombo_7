# Microsoft Developer Studio Project File - Name="main" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=main - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "main.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "main.mak" CFG="main - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "main - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "main - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName "Code"
# PROP Scc_LocalPath "."
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "main - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ""
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN" /D "NLIST" /D "DLIST" /D "COULOMB" /D "SHIFTF" /D "STYPE" /D "S14" /D "TRR" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib Ws2_32.lib /nologo /subsystem:console /profile /machine:I386

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN" /D "NLIST" /D "DLIST" /D "STYPE" /D "COULOMB" /D "TRR" /D "GOLIK" /D "DNA_GOLIK" /D "GOBT" /D "FSH" /D "TLATE" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib Ws2_32.lib /nologo /subsystem:console /profile /debug /machine:I386

!ENDIF 

# Begin Target

# Name "main - Win32 Release"
# Name "main - Win32 Debug"
# Begin Group "Read"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\allocate.cpp
# End Source File
# Begin Source File

SOURCE=.\read_atom.cpp
# End Source File
# Begin Source File

SOURCE=.\read_config.cpp
# End Source File
# Begin Source File

SOURCE=.\read_cons.cpp
# End Source File
# Begin Source File

SOURCE=.\read_dnagolik.cpp
# End Source File
# Begin Source File

SOURCE=.\read_golik.cpp
# End Source File
# Begin Source File

SOURCE=.\read_rest.cpp
# End Source File
# Begin Source File

SOURCE=.\read_sim.cpp
# End Source File
# Begin Source File

SOURCE=.\read_swap.cpp
# End Source File
# Begin Source File

SOURCE=.\read_topology.cpp
# End Source File
# Begin Source File

SOURCE=.\read_veloc.cpp
# End Source File
# Begin Source File

SOURCE=.\read_wall.cpp
# End Source File
# End Group
# Begin Group "Integrator"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\integrate.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_dos_mts.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_mts.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_mts_nhc.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_mts_npt_full.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_nhc.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_nhcp_full.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_npt_full.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_vverlet.cpp
# End Source File
# End Group
# Begin Group "Energy"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\calcvalue.cpp
# End Source File
# Begin Source File

SOURCE=.\cbend.cpp
# End Source File
# Begin Source File

SOURCE=.\cbond.cpp
# End Source File
# Begin Source File

SOURCE=.\cgolik.cpp
# End Source File
# Begin Source File

SOURCE=.\cimproper.cpp
# End Source File
# Begin Source File

SOURCE=.\cnonbond_nblist.cpp
# End Source File
# Begin Source File

SOURCE=.\crestraint.cpp
# End Source File
# Begin Source File

SOURCE=.\csasa.cpp
# End Source File
# Begin Source File

SOURCE=.\ctorsion.cpp
# End Source File
# Begin Source File

SOURCE=.\curey.cpp
# End Source File
# Begin Source File

SOURCE=.\cwall.cpp
# End Source File
# Begin Source File

SOURCE=.\ewald.cpp
# End Source File
# Begin Source File

SOURCE=.\kinet.cpp
# End Source File
# End Group
# Begin Group "Forces"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\force_long.cpp
# End Source File
# Begin Source File

SOURCE=.\force_short.cpp
# End Source File
# Begin Source File

SOURCE=.\forces.cpp
# End Source File
# End Group
# Begin Group "Output"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\movie.cpp
# End Source File
# Begin Source File

SOURCE=.\ofile.cpp
# End Source File
# Begin Source File

SOURCE=.\outend.cpp
# End Source File
# Begin Source File

SOURCE=.\output.cpp
# End Source File
# Begin Source File

SOURCE=.\save_config.cpp
# End Source File
# Begin Source File

SOURCE=.\save_veloc.cpp
# End Source File
# Begin Source File

SOURCE=.\xyz_config.cpp
# End Source File
# End Group
# Begin Group "Simulation"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\hybrid_MC.cpp
# End Source File
# Begin Source File

SOURCE=.\npt_md.cpp
# End Source File
# Begin Source File

SOURCE=.\nve_md.cpp
# End Source File
# Begin Source File

SOURCE=.\nvt_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\nvt_md.cpp
# End Source File
# Begin Source File

SOURCE=.\rem_opt.cpp
# End Source File
# Begin Source File

SOURCE=.\repexch.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_box.cpp
# End Source File
# End Group
# Begin Group "Tools"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\boxscale.cpp
# End Source File
# Begin Source File

SOURCE=.\curr_status.cpp
# End Source File
# Begin Source File

SOURCE=.\drift.cpp
# End Source File
# Begin Source File

SOURCE=.\end_to_end.cpp
# End Source File
# Begin Source File

SOURCE=.\helix.cpp
# End Source File
# Begin Source File

SOURCE=.\hessian.cpp
# End Source File
# Begin Source File

SOURCE=.\native_cont.cpp
# End Source File
# Begin Source File

SOURCE=.\quatfit.cpp
# End Source File
# Begin Source File

SOURCE=.\rad_gyr.cpp
# End Source File
# Begin Source File

SOURCE=.\rotate_xyz.cpp
# End Source File
# Begin Source File

SOURCE=.\translate.cpp
# End Source File
# Begin Source File

SOURCE=.\wall_angle.cpp
# End Source File
# Begin Source File

SOURCE=.\weight.cpp
# End Source File
# End Group
# Begin Group "DOS"

# PROP Default_Filter ""
# Begin Group "WLDOS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\dos.cpp
# End Source File
# Begin Source File

SOURCE=.\dos_assoc_trans.cpp
# End Source File
# Begin Source File

SOURCE=.\pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_dos.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_mpi.cpp
# End Source File
# Begin Source File

SOURCE=.\trans.cpp
# End Source File
# End Group
# Begin Group "MMDOS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\mmdos.cpp
# End Source File
# Begin Source File

SOURCE=.\mmdos_dke.cpp
# End Source File
# Begin Source File

SOURCE=.\mmdos_hmc.cpp
# End Source File
# Begin Source File

SOURCE=.\mmdos_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_mm_dos.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_mm_mpi.cpp
# End Source File
# End Group
# Begin Group "CTDOS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\ctdos.cpp
# End Source File
# Begin Source File

SOURCE=.\ctdos_assoc_trans.cpp
# End Source File
# Begin Source File

SOURCE=.\ctdos_hmc.cpp
# End Source File
# Begin Source File

SOURCE=.\ctdos_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\ctdos_trans.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_ct_dos.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_ct_mpi.cpp
# End Source File
# End Group
# Begin Group "XEDOS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\avrg_xe_mpi.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_xe_dos.cpp
# End Source File
# Begin Source File

SOURCE=.\swap_xe_mpi.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_axial.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_cosol.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_hmc.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_linear.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_scale.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_solv.cpp
# End Source File
# Begin Source File

SOURCE=.\xedos_trans.cpp
# End Source File
# End Group
# Begin Group "FXEDOS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\fx_edos.cpp
# End Source File
# Begin Source File

SOURCE=.\fx_edos_hmc.cpp
# End Source File
# Begin Source File

SOURCE=.\fx_edos_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\fx_edos_trans.cpp
# End Source File
# Begin Source File

SOURCE=.\update_weights.cpp
# End Source File
# End Group
# Begin Group "TDXEDOS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\avrg_tdxe_mpi.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_axial.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_hmc.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_init.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_linear.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_scale.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_solv.cpp
# End Source File
# Begin Source File

SOURCE=.\tdxedos_trans.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=.\dos_awrite.cpp
# End Source File
# Begin Source File

SOURCE=.\dos_interpolate.cpp
# End Source File
# Begin Source File

SOURCE=.\dos_stats.cpp
# End Source File
# Begin Source File

SOURCE=.\dos_svalues.cpp
# End Source File
# Begin Source File

SOURCE=.\dos_vblock.cpp
# End Source File
# Begin Source File

SOURCE=.\dos_wwrite.cpp
# End Source File
# End Group
# Begin Group "Minimization"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\brent.cpp
# End Source File
# Begin Source File

SOURCE=.\f1dim.cpp
# End Source File
# Begin Source File

SOURCE=.\frprmn.cpp
# End Source File
# Begin Source File

SOURCE=.\linmin.cpp
# End Source File
# Begin Source File

SOURCE=.\min_nma.cpp
# End Source File
# Begin Source File

SOURCE=.\minim.cpp
# End Source File
# Begin Source File

SOURCE=.\mnbrak.cpp
# End Source File
# Begin Source File

SOURCE=.\nrutil.cpp
# End Source File
# End Group
# Begin Group "Utilities"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\checks.cpp
# End Source File
# Begin Source File

SOURCE=.\clist.cpp
# End Source File
# Begin Source File

SOURCE=.\dgauss.cpp
# End Source File
# Begin Source File

SOURCE=.\min_image.cpp
# End Source File
# Begin Source File

SOURCE=.\nblist.cpp
# End Source File
# Begin Source File

SOURCE=.\nblist_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\pbc_all.cpp
# End Source File
# Begin Source File

SOURCE=.\pbc_chk.cpp
# End Source File
# Begin Source File

SOURCE=.\pbc_init.cpp
# End Source File
# Begin Source File

SOURCE=.\pbc_npt_full.cpp
# End Source File
# Begin Source File

SOURCE=.\ran_int.cpp
# End Source File
# Begin Source File

SOURCE=.\random.cpp
# End Source File
# Begin Source File

SOURCE=.\tak_histogram.cpp
# End Source File
# Begin Source File

SOURCE=.\tak_histogram.h
# End Source File
# Begin Source File

SOURCE=.\twham.cpp
# End Source File
# Begin Source File

SOURCE=.\twham.h
# End Source File
# Begin Source File

SOURCE=.\velscale.cpp
# End Source File
# Begin Source File

SOURCE=.\vinit.cpp
# End Source File
# Begin Source File

SOURCE=.\vinit_hyb.cpp
# End Source File
# End Group
# Begin Group "ConfigT"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\cbend2.cpp
# End Source File
# Begin Source File

SOURCE=.\cbond2.cpp
# End Source File
# Begin Source File

SOURCE=.\cgolik2.cpp
# End Source File
# Begin Source File

SOURCE=.\cimproper2.cpp
# End Source File
# Begin Source File

SOURCE=.\cnonbond_nblist2.cpp
# End Source File
# Begin Source File

SOURCE=.\configtemp.cpp
# End Source File
# Begin Source File

SOURCE=.\crestraint2.cpp
# End Source File
# Begin Source File

SOURCE=.\csasa2.cpp
# End Source File
# Begin Source File

SOURCE=.\ctorsion2.cpp
# End Source File
# Begin Source File

SOURCE=.\curey2.cpp
# End Source File
# Begin Source File

SOURCE=.\cwall2.cpp
# End Source File
# Begin Source File

SOURCE=.\ewald2.cpp
# End Source File
# End Group
# Begin Group "SMD"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\csmd.cpp
# End Source File
# Begin Source File

SOURCE=.\integrate_smd.cpp
# End Source File
# Begin Source File

SOURCE=.\smd_work.cpp
# End Source File
# End Group
# Begin Group "Misc"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\boxinv.cpp
# End Source File
# Begin Source File

SOURCE=.\cplx.cpp
# End Source File
# Begin Source File

SOURCE=.\init.cpp
# End Source File
# Begin Source File

SOURCE=.\isokin.cpp
# End Source File
# Begin Source File

SOURCE=.\ljpot.cpp
# End Source File
# Begin Source File

SOURCE=.\matrixmath.cpp
# End Source File
# Begin Source File

SOURCE=.\setlj.cpp
# End Source File
# Begin Source File

SOURCE=.\svalues.cpp
# End Source File
# Begin Source File

SOURCE=.\vblock.cpp
# End Source File
# Begin Source File

SOURCE=.\virialcor.cpp
# End Source File
# End Group
# Begin Group "MonteCarlo"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\ebend_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\ebond_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\egolik_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\eimpr_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\enbnd_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\erest_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\esasa_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\etors_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\eurey_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\ewald_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\forces_mc.cpp
# End Source File
# End Group
# Begin Group "Gromacs"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\write_gro.cpp
# End Source File
# Begin Source File

SOURCE=.\write_param.cpp
# End Source File
# Begin Source File

SOURCE=.\write_trr.cpp
# End Source File
# End Group
# Begin Group "MCMoves"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\mcmove_bkp.cpp
# End Source File
# Begin Source File

SOURCE=.\mcmove_displ.cpp
# End Source File
# Begin Source File

SOURCE=.\mcmove_hmc.cpp
# End Source File
# Begin Source File

SOURCE=.\mcmove_pivot.cpp
# End Source File
# Begin Source File

SOURCE=.\mcmove_trans.cpp
# End Source File
# End Group
# Begin Group "MonteCarlo2"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\ebend2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\ebond2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\egolik2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\eimpr2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\enbnd2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\erest2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\esasa2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\etors2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\eurey2_mc.cpp
# End Source File
# Begin Source File

SOURCE=.\ewald2_mc.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=.\defines.h
# End Source File
# Begin Source File

SOURCE=.\ioflush.cpp
# End Source File
# Begin Source File

SOURCE=.\main.cpp

!IF  "$(CFG)" == "main - Win32 Release"

# ADD CPP /Ob1
# SUBTRACT CPP /u

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\nr.h
# End Source File
# Begin Source File

SOURCE=.\nrutil.h
# End Source File
# End Target
# End Project
