import FWCore.ParameterSet.Config as cms

process = cms.Process("DDCMSDetectorTest")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.DDDetectorESProducer = cms.ESSource("DDDetectorESProducer",
                                            confGeomXMLFiles = cms.FileInPath('DetectorDescription/DDCMS/data/cms-tracker.xml'),
                                            appendToDataLabel = cms.string('CMS')
                                            )

process.DDVectorRegistryESProducer = cms.ESProducer("DDVectorRegistryESProducer",
                                                    appendToDataLabel = cms.string('CMS')
                                                    )

process.test = cms.EDAnalyzer("DDCMSDetector",
                              DDDetector = cms.ESInputTag('CMS')
                              )

process.testVectors = cms.EDAnalyzer("DDTestVectors",
                                     DDDetector = cms.ESInputTag('CMS')
                                     )

process.testDump = cms.EDAnalyzer("DDTestDumpFile",
                                  DDDetector = cms.ESInputTag('CMS')
                                  )

process.p = cms.Path(process.test+process.testVectors+process.testDump)
