      integer num1Dbins
      parameter(num1Dbins=n_p_el*n_th_el)
      real fprobsumm(num1Dbins)	
      integer thBins(num1Dbins),pBins(num1Dbins)
      common/csbnum/thBins,pBins
      common/cstsum/fprobsumm
