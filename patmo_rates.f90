module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,ntot(cellsNumber)
    integer::icell

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !O + O2 -> O3
      krate(icell,1) = 3.11d-34*(T/298)**(-2.0)*ntot(icell)

      !O + O3 -> O2 + O2
      krate(icell,2) = 1.83d-11*exp(-2164/T)

      !H + HO2 -> O + H2O
      krate(icell,3) = 1.6d-12

      !O + CO -> CO2
      krate(icell,4) = 1.7d-33*exp(-1510/T)*ntot(icell)

      !H + O2 -> HO2
      krate(icell,5) = (4.4d-32*(T/300)**(-1.3d0))*ntot(icell)/(1d0+((4.4d-32*(T/300)**(-1.3d0))*ntot(icell)/(7.5d-11*(T/300)**(2d-1))))*6d-1**((1d0+(log10((4.4d-32*(T/300)**(-1.3d0))*ntot(icell)/(7.5d-11*(T/300)**(2d-1))))**(2d0))**(-1d0))

    end do

  end subroutine computeRates

end module patmo_rates
