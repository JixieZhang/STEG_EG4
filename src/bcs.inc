c      include file BCS
c
c_begin_doc
c  RCS ID string
c  $Id: bcs.inc 1.6 1996/01/05 22:07:41 freyberg Exp $
c
c  Documentation for include file BCS
c
c  Purpose: Standard common for BOS banks
c  --------
c
c  Used in routines: all
c  ----------------
c
c  Notes:
c  ------
c
c  Author:   Dieter Cords      Created:  Mon Jun 12 17:50:51 EDT 1995
c  -------
c
c  RCS information: 
      CHARACTER*(132)  CIBCS
      DATA CIBCS/
     1'$Id: bcs.inc 1.6 1996/01/05 22:07:41 freyberg Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in BCS common block
c  --------------------------------------------
c
c_end_var
c
      INTEGER Nbcs
      PARAMETER (Nbcs = 700000)
c
      INTEGER IWW(Nbcs), JUNK(5)
      REAL RW(Nbcs)
      INTEGER*2 IW16(2*Nbcs)
      EQUIVALENCE (IWW(1),RW(1),IW16(1))
c 
c 
c IWCOL and IWROW provide a quick way to find the number of columns
c and rows for a particular bank
c
c
      COMMON/BCS/JUNK,IWW
      INTEGER IWROW(Nbcs), IWCOL(Nbcs), IWNAME(Nbcs), IWNR(Nbcs), 
     1        IWNXT(Nbcs)
      EQUIVALENCE (IWW(1), IWCOL(6))
      EQUIVALENCE (IWW(1), IWROW(5))
      EQUIVALENCE (IWW(1), IWNAME(4))
      EQUIVALENCE (IWW(1), IWNR(3))
      EQUIVALENCE (IWW(1), IWNXT(2))
c
c
      save /BCS/
c----6----------------------------------------------------------------72

