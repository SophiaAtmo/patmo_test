module patmo_commons
  implicit none

  integer,parameter::reactionsNumber = 13
  integer,parameter::chemReactionsNumber = 5
  integer,parameter::photoReactionsNumber = 3
  integer,parameter::reverseReactionsNumber = 5
  integer,parameter::chemSpeciesNumber = 9
  integer,parameter::speciesNumber = 11
  integer,parameter::positionTgas = 10
  integer,parameter::positionDummy = 11
  integer,parameter::cellsNumber = 100
  integer,parameter::photoBinsNumber = 6800
  integer,parameter::patmo_idx_O = 1
  integer,parameter::patmo_idx_O2 = 2
  integer,parameter::patmo_idx_M = 3
  integer,parameter::patmo_idx_O3 = 4
  integer,parameter::patmo_idx_H = 5
  integer,parameter::patmo_idx_HO2 = 6
  integer,parameter::patmo_idx_H2O = 7
  integer,parameter::patmo_idx_CO = 8
  integer,parameter::patmo_idx_CO2 = 9

  integer,parameter::chemReactionsOffset = 0
  integer,parameter::photoReactionsOffset = chemReactionsNumber
  integer,parameter::reverseReactionsOffset = &
      photoReactionsOffset + photoReactionsNumber

  integer,parameter::neqAll = speciesNumber*cellsNumber
  integer,parameter::maxNameLength = 50

  integer,dimension(photoReactionsNumber)::photoPartnerIndex = (/patmo_idx_O2,patmo_idx_O3,patmo_idx_CO2/)

  integer,parameter,dimension(reactionsNumber)::indexReactants1 = (/patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_CO2,&
      patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_CO2,&
      patmo_idx_HO2/)
  integer,parameter,dimension(reactionsNumber)::indexReactants2 = (/patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_HO2,&
      patmo_idx_CO,&
      patmo_idx_O2,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      patmo_idx_M,&
      patmo_idx_O2,&
      patmo_idx_H2O,&
      patmo_idx_M,&
      patmo_idx_M/)
  integer,parameter,dimension(reactionsNumber)::indexReactants3 = (/patmo_idx_M,&
      positionDummy,&
      positionDummy,&
      patmo_idx_M,&
      patmo_idx_M,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy,&
      positionDummy/)
  integer,parameter,dimension(reactionsNumber)::indexProducts1 = (/patmo_idx_O3,&
      patmo_idx_O2,&
      patmo_idx_O,&
      patmo_idx_CO2,&
      patmo_idx_HO2,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_CO,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_H,&
      patmo_idx_O,&
      patmo_idx_H/)
  integer,parameter,dimension(reactionsNumber)::indexProducts2 = (/patmo_idx_M,&
      patmo_idx_O2,&
      patmo_idx_H2O,&
      patmo_idx_M,&
      patmo_idx_M,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_HO2,&
      patmo_idx_CO,&
      patmo_idx_O2/)

end module patmo_commons
