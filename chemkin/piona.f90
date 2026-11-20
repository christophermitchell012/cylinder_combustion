! Compare a binary C5 vapor using the legacy double-precision CHEMKIN-II API.
! Read kelvin values from stdin, one per line. CSV goes to stdout.
program piona
  use iso_fortran_env, only: real64, input_unit, output_unit, error_unit, iostat_end
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: nInt, nReal, nChar, linkUnit, status, nElem, nSpec, nRxn, nFit
  integer, allocatable :: iWork(:)
  real(real64), allocatable :: rWork(:)
  character(16), allocatable :: cWork(:)
  character(16) :: names(2)
  logical :: badNames
  real(real64) :: tempK, gMol(2), gasR, gasRcal, atm, eqIso, outletIso, shift
  real(real64), parameter :: feedIso=0.20_real64, approach=0.80_real64, blend=0.10_real64

  open(newunit=linkUnit, file='chem.bin', status='old', action='read', form='unformatted')
  call cklen(linkUnit,error_unit,nInt,nReal,nChar)
  allocate(iWork(nInt),rWork(nReal),cWork(nChar))
  call ckinit(nInt,nReal,nChar,linkUnit,error_unit,iWork,rWork,cWork)
  close(linkUnit)
  call ckindx(iWork,rWork,nElem,nSpec,nRxn,nFit)
  if (nSpec /= 2) error stop 'Expected the supplied two-species deck'
  call cksyms(cWork,error_unit,names,badNames)
  if (badNames .or. any(names /= ['NC5H12','IC5H12'])) error stop 'Check species order'
  ! CKGML uses ergs/mol. Use the library gas constant in matching units.
  call ckrp(iWork,rWork,gasR,gasRcal,atm)
  write(output_unit,'(a)') 'temperature_K,equilibrium_iso,outlet_iso,delta_P_pp,delta_I_pp'
  do
    read(input_unit,*,iostat=status) tempK
    if (status == iostat_end) exit
    if (status /= 0) error stop 'Expected a numeric temperature in kelvin'
    if (.not.ieee_is_finite(tempK) .or. tempK < 298.15_real64 .or. tempK > 5000) &
        error stop 'Outside common NASA7 range: 298.15 to 5000 K'
    call ckgml(tempK,iWork,rWork,gMol)
    ! Both species have the same molecular mass, so this is also a mass fraction.
    eqIso = 1 / (1 + exp((gMol(2)-gMol(1))/(gasR*tempK)))
    outletIso = feedIso + approach*(eqIso-feedIso)
    shift = 100*blend*(outletIso-feedIso)
    write(output_unit,'(es24.16,4(",",es24.16))') tempK,eqIso,outletIso,-shift,shift
  end do
end program
