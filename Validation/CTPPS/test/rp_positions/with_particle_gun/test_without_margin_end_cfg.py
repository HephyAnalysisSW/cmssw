import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process('CTPPSFastSimulation', eras.ctpps_2016)

# minimal logger settings
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# particle-data table
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# particle generator
process.generator = cms.EDProducer("RandomtXiGunProducer",
  Verbosity = cms.untracked.int32(0),

  FireBackward = cms.bool(True),
  FireForward = cms.bool(True),

  PGunParameters = cms.PSet(
    PartID = cms.vint32(2212),
    ECMS = cms.double(13E3),

    Mint = cms.double(0),
    Maxt = cms.double(2),
    MinXi = cms.double(0.0),
    MaxXi = cms.double(0.2),

    MinPhi = cms.double(-3.14159265359),
    MaxPhi = cms.double(+3.14159265359)
  )
)

# random seeds
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    sourceSeed = cms.PSet(initialSeed =cms.untracked.uint32(98765)),
    generator = cms.PSet(initialSeed = cms.untracked.uint32(98766)),
    SmearingGenerator = cms.PSet(initialSeed =cms.untracked.uint32(3849))
)

# geometry
from Geometry.VeryForwardGeometry.geometryRP_cfi import totemGeomXMLFiles, ctppsDiamondGeomXMLFiles

process.XMLIdealGeometryESSource_CTPPS = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = totemGeomXMLFiles+ctppsDiamondGeomXMLFiles+["Validation/CTPPS/test/rp_positions/data/2016_preTS2_without_margin_end/RP_Dist_Beam_Cent.xml"],
    rootNodeName = cms.string('cms:CMSE'),
)

process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

# fast simulation
process.load('SimCTPPS.OpticsParameterisation.ctppsFastProtonSimulation_cfi')
process.ctppsFastProtonSimulation.produceScoringPlaneHits = False
process.ctppsFastProtonSimulation.produceRecHits = True
process.ctppsFastProtonSimulation.checkApertures = False
process.ctppsFastProtonSimulation.produceHitsRelativeToBeam = False
process.ctppsFastProtonSimulation.roundToPitch = True

# strips reco: pattern recognition
process.load('RecoCTPPS.TotemRPLocal.totemRPUVPatternFinder_cfi')
process.totemRPUVPatternFinder.tagRecHit = cms.InputTag('ctppsFastProtonSimulation')
process.ctppsFastProtonSimulation.produceHitsRelativeToBeam = True

# strips reco: track fitting
process.load('RecoCTPPS.TotemRPLocal.totemRPLocalTrackFitter_cfi')

# common reco: lite track production
process.load('RecoCTPPS.TotemRPLocal.ctppsLocalTrackLiteProducer_cfi')

# distribution plotter
process.ctppsTrackDistributionPlotter_reco = cms.EDAnalyzer("CTPPSTrackDistributionPlotter",
  tracksTag = cms.InputTag("ctppsLocalTrackLiteProducer"),
  outputFile = cms.string("output_without_margin_end.root")
)

# processing path
process.p = cms.Path(
    process.generator

    * process.ctppsFastProtonSimulation

    * process.totemRPUVPatternFinder
    * process.totemRPLocalTrackFitter
    * process.ctppsLocalTrackLiteProducer

    * process.ctppsTrackDistributionPlotter_reco
)
