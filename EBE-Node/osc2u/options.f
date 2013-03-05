c $Id: options.f,v 1.5 1997/08/25 08:17:21 weber Exp $
c... law: include file (only) for global parameters & options
      integer numcto,numctp
      parameter(numcto=400,numctp=400)
c...
      integer   CTOption(numcto)
      character ctodc(numcto)*2
c...
      real*8    CTParam(numctp)
      character ctpdc(numctp)*2

      logical bf13,bf14,bf15,bf16,bf17,bf18,bf19,fixedseed
      common /options/CTOption,CTParam
      common /optstrings/ctodc,ctpdc
      common /loptions/fixedseed,bf13,bf14,bf15,bf16,bf17,bf18,
     .     bf19
