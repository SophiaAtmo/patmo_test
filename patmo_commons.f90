module patmo_commons
  implicit none

  integer,parameter::reactionsNumber = 13
  integer,parameter::chemReactionsNumber = 5
  integer,parameter::photoReactionsNumber = 3
  integer,parameter::reverseReactionsNumber = 5
  integer,parameter::chemSpeciesNumber = 8
  integer,parameter::speciesNumber = 10
  integer,parameter::positionTgas = 9
  integer,parameter::positionDummy = 10
  integer,parameter::cellsNumber = 100
  integer,parameter::photoBinsNumber = 6800
  integer,parameter::patmo_idx_O = 1
  integer,parameter::patmo_idx_O2 = 2
  integer,parameter::patmo_idx_O3 = 3
  integer,parameter::patmo_idx_H = 4
  integer,parameter::patmo_idx_HO2 = 5
  integer,parameter::patmo_idx_H2O = 6
  integer,parameter::patmo_idx_CO = 7
  integer,parameter::patmo_idx_CO2 = 8

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
      positionDummy,&
      patmo_idx_O2,&
      patmo_idx_H2O,&
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
  integer,parameter,dimension(reactionsNumber)::indexProducts2 = (/positionDummy,&
      patmo_idx_O2,&
      patmo_idx_H2O,&
      positionDummy,&
      positionDummy,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O,&
      patmo_idx_O2,&
      patmo_idx_O3,&
      patmo_idx_HO2,&
      patmo_idx_CO,&
      patmo_idx_O2/)

end module patmo_commons
