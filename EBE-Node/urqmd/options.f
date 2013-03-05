c $Id: options.f,v 1.6 1998/06/15 13:35:28 weber Exp $
c... law: include file (only) for global parameters & options
      integer numcto,numctp
      parameter(numcto=400) ! maximum number of options
      parameter(numctp=400) ! maximum number of parameters
c...
      integer   CTOption(numcto)
      character ctodc(numcto)*2
c...
      real*8    CTParam(numctp)
      character ctpdc(numctp)*2

      logical bf13,bf14,bf15,bf16,bf17,bf18,bf19,fixedseed
      common /options/CTOption,CTParam
      common /optstrings/ctodc,ctpdc
      common /loptions/fixedseed,bf13,bf14,bf15,bf16,bf17,bf18,bf19
