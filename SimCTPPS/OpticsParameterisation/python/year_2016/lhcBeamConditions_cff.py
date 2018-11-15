import FWCore.ParameterSet.Config as cms

lhcBeamConditions = cms.PSet(
    sqrtS = cms.double(13.e3), # in GeV
    vertexSize = cms.double(10.e-6), # in m
    beamDivergence = cms.double(20.e-6), # in rad

    # vertex offset in both sectors
    xOffsetSector45 = cms.double(0), # in m
    xOffsetSector56 = cms.double(0), # in m
    yOffsetSector45 = cms.double(300e-6), # in m
    yOffsetSector56 = cms.double(200e-6), # in m

    # crossing angle
    halfCrossingAngleSector45 = cms.double(179.394e-6), # in rad
    halfCrossingAngleSector56 = cms.double(191.541e-6), # in rad
)
