import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')
options.register('outfile',  '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, '') 

options.parseArguments()

print "Writing CTPPS EDM output to", options.outfile

from Configuration.StandardSequences.Eras import eras
process = cms.Process('CTPPSFastSimulation', eras.ctpps_2016)

#-- Message Logger ------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
  wantSummary = cms.untracked.bool(False),
  allowUnscheduled = cms.untracked.bool( True )
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.destinations = cms.untracked.vstring('cerr')
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(100)
)
process.MessageLogger.cout  = process.MessageLogger.cerr

# number of events
#process.source = cms.Source("EmptySource")
process.source = cms.Source('PoolSource',
                            noEventSort = cms.untracked.bool(True),                                        # add this                                                                       
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'), # and this    
                            fileNames = cms.untracked.vstring( options.inputFiles ) )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# particle-data table
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# particle generator
process.generator = cms.EDProducer("RandomXiThetaGunProducer",
  particleId = cms.uint32(2212),

  energy = cms.double(6500),  # nominal beam energy, GeV
  xi_min = cms.double(0.),
  xi_max = cms.double(0.25),
  theta_x_mean = cms.double(0),
  theta_x_sigma = cms.double(50E-6), # in rad
  theta_y_mean = cms.double(0),
  theta_y_sigma = cms.double(50E-6),

  nParticlesSector45 = cms.uint32(1),
  nParticlesSector56 = cms.uint32(1),
)

# random seeds
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    sourceSeed = cms.PSet(initialSeed =cms.untracked.uint32(98765)),
    generator = cms.PSet(initialSeed = cms.untracked.uint32(98766)),
    SmearingGenerator = cms.PSet(initialSeed =cms.untracked.uint32(3849)),
    beamDivergenceVtxGenerator = cms.PSet(initialSeed =cms.untracked.uint32(3849))  
)

# geometry
process.load("Geometry.VeryForwardGeometry.geometryRP_2017_cfi")
del(process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles[-1])
process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.append("Validation/CTPPS/test_2017/RP_Dist_Beam_Cent.xml")

# fast simulation
process.load('SimCTPPS.OpticsParameterisation.year_2017_OF.ctppsFastProtonSimulation_cfi')
#process.ctppsFastProtonSimulation.hepMCTag = cms.InputTag('generatorSmeared')
process.ctppsFastProtonSimulation.hepMCTag = cms.InputTag('beamDivergenceVtxGenerator')
#process.ctppsFastProtonSimulation.hepMCTag = cms.InputTag('source')
process.ctppsFastProtonSimulation.xangle = 0
process.ctppsFastProtonSimulation.produceScoringPlaneHits = False
process.ctppsFastProtonSimulation.produceRecHits = True
process.ctppsFastProtonSimulation.checkApertures = False
process.ctppsFastProtonSimulation.useEmpiricalApertures = True 
process.ctppsFastProtonSimulation.produceHitsRelativeToBeam = True 
process.ctppsFastProtonSimulation.roundToPitch = True

# beam-smearing settings
process.load("IOMC.EventVertexGenerators.beamDivergenceVtxGenerator_cfi")
#process.beamDivergenceVtxGenerator.src = cms.InputTag("generator", "unsmeared")
#process.beamDivergenceVtxGenerator.src = cms.InputTag("source", "")
process.beamDivergenceVtxGenerator.src = cms.InputTag('generatorSmeared')

process.beamDivergenceVtxGenerator.simulateBeamDivergence = True
process.beamDivergenceVtxGenerator.simulateVertex = True

# values in rad
process.beamDivergenceVtxGenerator.beamDivergenceX = 20E-6
process.beamDivergenceVtxGenerator.beamDivergenceY = 20E-6

# values in cm
process.beamDivergenceVtxGenerator.vertexMeanX = 0.02476
process.beamDivergenceVtxGenerator.vertexMeanY = -0.06920
process.beamDivergenceVtxGenerator.vertexMeanZ = -0.8775

process.beamDivergenceVtxGenerator.vertexSigmaX = 0.
process.beamDivergenceVtxGenerator.vertexSigmaY = 0.
process.beamDivergenceVtxGenerator.vertexSigmaZ = 0.



# local track reco
process.load('RecoCTPPS.TotemRPLocal.totemRPUVPatternFinder_cfi')
process.totemRPUVPatternFinder.tagRecHit = cms.InputTag('ctppsFastProtonSimulation')

process.load('RecoCTPPS.TotemRPLocal.totemRPLocalTrackFitter_cfi')

process.load("RecoCTPPS.PixelLocal.ctppsPixelLocalTracks_cfi")
process.ctppsPixelLocalTracks.label = "ctppsFastProtonSimulation"

process.load('RecoCTPPS.TotemRPLocal.ctppsLocalTrackLiteProducer_cff')
process.ctppsLocalTrackLiteProducer.includeDiamonds = False


process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string(options.outfile),
    outputCommands = cms.untracked.vstring([ 'drop *', 'keep *_*_*_CTPPSFastSimulation', 'keep *_slimmedJets_*_*' ]) 
)

process.load("RecoCTPPS.ProtonReconstruction.year_2017_OF.ctppsProtonReconstructionOF_cfi")
process.ctppsProtonReconstructionOFDB.applyExperimentalAlignment = False # do not use alignment for LHC data

process.ctppsLHCInfoESSource = cms.ESSource("CTPPSLHCInfoESSource",
  beamEnergy = cms.double(6500),
  xangle = cms.double(0)
)

process.reco_step = cms.Path(
  process.ctppsProtonReconstructionOFDB ##<- need to run on all data & MC in mAOD
  )

process.dump=cms.EDAnalyzer('EventContentAnalyzer')
# processing path
process.p = cms.Path(
  #process.dump*
    process.beamDivergenceVtxGenerator*
  process.ctppsFastProtonSimulation
)
process.simulation_step = cms.Path(
     process.totemRPUVPatternFinder
    * process.totemRPLocalTrackFitter
    * process.ctppsPixelLocalTracks
    * process.ctppsLocalTrackLiteProducer
)
process.outpath = cms.EndPath(process.out)

process.schedule = cms.Schedule(
#    process.p,
#    process.simulation_step,
    process.reco_step,
    process.outpath
)


def UseCrossingAngle150():
  process.ctppsFastProtonSimulation.xangle = process.ctppsLHCInfoESSource.xangle = 150
  process.ctppsFastProtonSimulation.empiricalAperture45_xi0 = 0.158
  process.ctppsFastProtonSimulation.empiricalAperture56_xi0 = 0.20
  #process.ctppsAcceptancePlotter.outputFile = "acceptance_xangle_150.root"

def UseCrossingAngle140():
  process.ctppsFastProtonSimulation.xangle = process.ctppsLHCInfoESSource.xangle = 140
  process.ctppsFastProtonSimulation.empiricalAperture45_xi0 = 0.153
  process.ctppsFastProtonSimulation.empiricalAperture56_xi0 = 0.19
  #process.ctppsAcceptancePlotter.outputFile = "acceptance_xangle_140.root"

def UseCrossingAngle130():
  process.ctppsFastProtonSimulation.xangle = process.ctppsLHCInfoESSource.xangle = 130
  process.ctppsFastProtonSimulation.empiricalAperture45_xi0 = 0.148
  process.ctppsFastProtonSimulation.empiricalAperture56_xi0 = 0.18
  #process.ctppsAcceptancePlotter.outputFile = "acceptance_xangle_130.root"

UseCrossingAngle130()
