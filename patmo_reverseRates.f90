module patmo_reverseRates
contains

  !compute reverse rates using thermochemistry polynomials
  subroutine computeReverseRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(:)
    real*8::Tgas(cellsNumber)
    real*8::lnTgas(cellsNumber)
    real*8::Tgas2(cellsNumber)
    real*8::Tgas3(cellsNumber)
    real*8::Tgas4(cellsNumber)
    real*8::invTgas(cellsNumber)
    real*8::ntot(cellsNumber)
    integer::i

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

    !extrapolate lower and upper limits
    do i=1,cellsNumber
      Tgas(i) = max(inTgas(i),2d2)
      Tgas(i) = min(Tgas(i),5d3)
    end do

    lnTgas(:) = log(Tgas(:))
    Tgas2(:) = Tgas(:)**2
    Tgas3(:) = Tgas(:)**3
    Tgas4(:) = Tgas(:)**4
    invTgas(:) = 1d0/Tgas(:)

    !O3 -> O + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,9) = krate(i,1)*exp(3.543341d0*(lnTgas(i)-1d0) &
            - 4.164922d-3*Tgas(i) &
            + 4.402935d-7*Tgas2(i) &
            + 5.434827d-10*Tgas3(i) &
            - 2.202172d-13*Tgas4(i) &
            - 1.219382d4*invTgas(i) &
            - 2.572867d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,9) = krate(i,1)*exp(-6.125694d0*(lnTgas(i)-1d0) &
            + 6.280764d-3*Tgas(i) &
            - 1.355459d-6*Tgas2(i) &
            + 1.4979d-10*Tgas3(i) &
            - 6.392726d-15*Tgas4(i) &
            - 1.533445d4*invTgas(i) &
            + 4.921999d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,9) = 0d0
      end if
      !divided because pseudo-3body
      krate(i,9) = krate(i,9) / ntot(i)
    end do

    !O2 + O2 -> O + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,10) = krate(i,2)*exp(-9.892634d-1*(lnTgas(i)-1d0) &
            + 2.38397d-3*Tgas(i) &
            + 1.328442d-7*Tgas2(i) &
            - 7.580525d-10*Tgas3(i) &
            + 2.692968d-13*Tgas4(i) &
            - 4.711464d4*invTgas(i) &
            + 3.019058d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,10) = krate(i,2)*exp(7.552007d0*(lnTgas(i)-1d0) &
            - 6.636263d-3*Tgas(i) &
            + 1.377587d-6*Tgas2(i) &
            - 1.506792d-10*Tgas3(i) &
            + 6.409727d-15*Tgas4(i) &
            - 4.433355d4*invTgas(i) &
            - 4.279077d1)
      else
        krate(i,10) = 0d0
      end if
    end do

    !O + H2O -> H + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,11) = krate(i,3)*exp(-5.651042d-1*(lnTgas(i)-1d0) &
            + 2.832998d-4*Tgas(i) &
            + 1.332481d-6*Tgas2(i) &
            - 1.055033d-9*Tgas3(i) &
            + 2.703812d-13*Tgas4(i) &
            - 2.690915d4*invTgas(i) &
            + 2.067055d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,11) = krate(i,3)*exp(1.451612d0*(lnTgas(i)-1d0) &
            - 5.323445d-4*Tgas(i) &
            + 7.194698d-8*Tgas2(i) &
            - 6.660215d-12*Tgas3(i) &
            + 2.462405d-16*Tgas4(i) &
            - 2.616456d4*invTgas(i) &
            - 9.293851d0)
      else
        krate(i,11) = 0d0
      end if
    end do

    !CO2 -> O + CO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,12) = krate(i,4)*exp(4.390988d0*(lnTgas(i)-1d0) &
            - 6.436901d-3*Tgas(i) &
            + 2.463657d-6*Tgas2(i) &
            - 6.398634d-10*Tgas3(i) &
            + 6.755604d-14*Tgas4(i) &
            - 6.315014d4*invTgas(i) &
            - 4.340561d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,12) = krate(i,4)*exp(9.556118d-1*(lnTgas(i)-1d0) &
            - 7.085225d-4*Tgas(i) &
            + 8.431887d-8*Tgas2(i) &
            - 6.381516d-12*Tgas3(i) &
            + 1.992179d-16*Tgas4(i) &
            - 6.39848d4*invTgas(i) &
            + 1.287429d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,12) = 0d0
      end if
      !divided because pseudo-3body
      krate(i,12) = krate(i,12) / ntot(i)
    end do

    !HO2 -> H + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,13) = krate(i,5)*exp(1.980658d0*(lnTgas(i)-1d0) &
            + 8.761934d-4*Tgas(i) &
            - 1.885165d-6*Tgas2(i) &
            + 1.216258d-9*Tgas3(i) &
            - 3.024262d-13*Tgas4(i) &
            - 2.41457d4*invTgas(i) &
            - 5.056693d-1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,13) = krate(i,5)*exp(1.988673d0*(lnTgas(i)-1d0) &
            - 6.124052d-4*Tgas(i) &
            + 3.418794d-8*Tgas2(i) &
            + 9.283655d-14*Tgas3(i) &
            - 7.376956d-17*Tgas4(i) &
            - 2.422666d4*invTgas(i) &
            + 1.100322d-2)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,13) = 0d0
      end if
      !divided because pseudo-3body
      krate(i,13) = krate(i,13) / ntot(i)
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates
