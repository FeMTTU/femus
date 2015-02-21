CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TimedepNavierStokes)

femusMacroBuildApplication(NavierStokes/BackwardEulerNavierStokesParallel backeuler)

