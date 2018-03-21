# STEG_EG4
By Jixie Zhang @20180306
This program is called STEG. It is used to generate xs map and generate events
for simulation. So far it is dedicated for EG4 experiment in Jlab.

---------------------------------------------------------------------
03/19/2018
1. Fixed typo in Makefile. very critical!
2. Changed stegAll_myOwnNoAvgXSinBinMultSc.f to make it compatible with Lamiaa's version.
   If file 'stegAllvPREG.dat' exists in currrent dir, will read arguments from this file.
   Otherwise, ir will read arguments from shell (usually run in this way: 
   stegAllvKPA_rhel7 < inputfile )
3. Changed stegAll_myOwnNoAvgXSinBinMultSc.f. During generating events, print event# 
   message every 1000 events. Before it prints this message every events.
4. Fixed bug in file nThnPnSum.inc. 
-      parameter(num1Dbins=n_p_el*n_p_el)
+      parameter(num1Dbins=n_p_el*n_th_el)
   This bug does not do anything bad in the past since n_th_el is equal to n_p_el. But it will
   if n_th_el is not equal to n_p_el.

---------------------------------------------------------------------
03/15/2018
1. Solved bug about signning EB_INDEX. 
2. Files newSFsSEK17.f sek_meKPA.f sigmarNEW_SEKnw17.f are in CRLF(windows) ending format, 
   Changed them to CR(Unix) ending format. Changed all of them to EMACS indent format.
3. Replaced YoniSEK.f using the original version without my format.
4. Added file Changes_Since_V1.0.txt make it possible to trace back to version 1.0. 

---------------------------------------------------------------------
1. Found that the xs under elastic peak seems abnormal large than its neighbour points.
This is due to the fact that in some theta bin, the center value of the theta-p bin is too close to the
E'_elastic. Since we are not simulate elastic peak in STEG, we would like to remove these points.
The following cuts are applied:

    XCUT < A*0.99; 
    E' < E'_elastic -E'_bin_width; 
    W > 0.95 .

2. During throw event, called Xiaochao's subroutine 'rleg4_simple(xxx)' to calculate TA,TB.
TA will be used to determine MSC for theta of the scattered electron. 

---------------------------------------------------------------------
03/06/2018
Reorganzied this program  and put it into github.

1. Essential files for successful compilation are all in src: *.f, *.inc
2. All essential file for running this program are in otherfiles:
 F2ALLM07.inp?, rcslacpol.file, *.ful, *.dat
3. Makefile changed   


---------------------------------------------------------------------
02/28/2018:
1. Added Xiaochao's subroutine 'rleg4_simple(xxx)' to calculate TA,TB,PF according to 
   given beam energy for NH3.
2. Changed set_things_up(xxx) to sign EB_INDEX according to given beam energy and target type.
3. Removed MSC in phi angle when throw events.
4. Changed XCUT to A*0.99D0. Previously it was 'XCUT = 1.0/(0.0658/(E*EP*sin2)+1.0)'.
5. Cleaned up some useless comments, removed some unused subroutine.

   TODO: clean up, reorganize, re-format then upload into github

---------------------------------------------------------------------
02/06/2018

1. STEG outputs 'infinite' or 'NAN' occasionally

According to fortran language manual, IEEE 754 floating point numbers can represent positive 
or negative infinity, and NaN (not a number). 
These three values arise from calculations whose result is undefined or cannot be represented 
accurately. They are defined as
  1/0 = ∞
  log (0) = -∞
  sqrt (-1) = NaN

I looked into the source code and the output. I found that these issues occur only at the edge of 
the kinematics acceptance, for example, E'>E_elas due to round off effect or due to not precise 
input value of Proton|Deuteron mass.  

	Feb 06: A subroutine to check infinity or NAN is added. if any of xsB, xsR, xsP and xsPR is 
	        infinite or NaN, set xsB=xsR=xsP=xsPR=0.0.
	Feb 06: solved!

2. STEG produces negative radiated unpolarized cross section, but born cross section are all 
   positive or zero. It happens only to NH3 target.

	Feb 06: Added code to check the radiated unpolarized xs, if negative, set xsB=xsR=xsP=xsPR=0.0.
	Feb 06: solved!

---------------------------------------------------------------------

10/13/2017: Debug for create xs map for NH3 target, use Sebastian's new code which can 
    alter FakeProton between TRUE and FALSE. 
	
Sep. 2017 Fixed bug in writing bos outpt

Sep. 2017 Installed CLAS package in CentOS7. Adapted both STEG and CombineFracMaps programs into my own environment.

June 2017 Copied Krishna's version 


---------------------------------------------------------------------
8/14/16:

Choose proper values of the four parameters/flags (changeA1, AsymChoiceKpa, SFchoiceKpa, and 
       factorRA) in kpVarChanges.inc depending on which simulation we want to perform.
Choose n_th_el=2, and n_p_el=200 etc in the file binning.inc to make exe for XS-map generation
Choose n_th_el=n_p_el=200 or something like that for exe to use for event sampling/generation

