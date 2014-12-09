############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")
 
 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')

## pixel additions
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
Source_Files = cms.untracked.vstring(
    ## single muons PU=0
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_1.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_2.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_3.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_4.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_5.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_6.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_7.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_8.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_9.root',
    '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Muons/NoPU/SingleMuon_DIGI_10.root',
    ## single muons PU=140
    #'/store/mc/TTI2023Upg14D/SingleMuMinusFlatPt0p2To150/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0097A61E-04E7-E311-824C-003048678F9C.root',
    )
process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleMuon_noPU_TrkPerf.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# Path definitions & schedule
############################################################

#run the tracking
BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)

#pixel stuff
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
process.siPixelRecHits = siPixelRecHits

process.L1PixelTrackFit = cms.EDProducer("L1PixelTrackFit")
process.pixTrk = cms.Path(process.L1PixelTrackFit)

process.pixRec = cms.Path(
    process.RawToDigi+
    process.siPixelRecHits
)
process.raw2digi_step = cms.Path(process.RawToDigi)


# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(13),
                                       DebugMode = cms.bool(False),      ## printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   ## save all L1 tracks, not just truth matched to primary particle
                                       DoPixelTrack = cms.bool(True)     ## save information for pixel tracks
                                       )
process.ana = cms.Path(process.L1TrackNtuple)

#output module (can use this to store edm-file instead of root-ntuple)
#process.out = cms.OutputModule( "PoolOutputModule",
#                                fileName = cms.untracked.string("FileOut.root"),
#                                fastCloning = cms.untracked.bool( False ),
#                                outputCommands = cms.untracked.vstring('drop *',
#                                                                       'keep *_*_Level1PixelTracks_*',
#                                                                       'keep *_*_Level1TTTracks_*',
#                                                                       'keep *_*_StubAccepted_*',
#                                                                       'keep *_*_ClusterAccepted_*',
#                                                                       'keep *_*_MergedTrackTruth_*')
#)
#process.FEVToutput_step = cms.EndPath(process.out)


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI
process = cust_2023TTI(process)

#process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.pixRec,process.pixTrk,process.FEVToutput_step)
process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.pixRec,process.pixTrk,process.ana)
