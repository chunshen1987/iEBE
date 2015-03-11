add_executable (superMC.e
  main.cpp 
  Bases.cpp 
  MCnucl.cpp 
  GlueDensity.cpp 
  MakeDensity.cpp
  KLNModel.cpp 
  OverLap.cpp 
  Largex.cpp 
  Regge96.cpp 
  rcBKfunc.cpp
  MathBasics.cpp 
  ParameterReader.cpp
  arsenal.cpp
  EOS.cpp
  GaussianNucleonsCal.cpp 
  NBD.cpp 
  RandomVariable.cpp
  TableFunction.cpp 
  Table.cpp
  HulthenFunc.cpp
  GaussianDistribution.cpp
  Particle.cpp
  )
target_link_libraries (superMC.e ${LIBS})

install(TARGETS superMC.e DESTINATION ${CMAKE_HOME_DIRECTORY})
