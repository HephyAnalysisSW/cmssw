''' FWLite example
'''
# Standard imports
import ROOT
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *
from math import sqrt, sin, cos

small = False

# example file
events = Events(['ctppsSim.root'])

#vector<CTPPSLocalTrackLite>           "ctppsLocalTrackLiteProducer"   ""                "RECO"
#edm::TriggerResults            "TriggerResults"            ""        "CTPPSFastSimulation"
#vector<reco::ProtonTrack>      "ctppsProtonReconstructionOFDB"   ""        "CTPPSFastSimulation"

# RECO
edmCollections = {
'ptracks':{'type':'vector<reco::ProtonTrack>', 'label': ( "ctppsProtonReconstructionOFDB" ) },
   }

# add handles
for k, v in edmCollections.iteritems():
    v['handle'] = Handle(v['type'])

nevents = 1 if small else events.size()

for i in range(nevents):
    events.to(i)

    eaux  = events.eventAuxiliary()

    # run/lumi/event
    run   = eaux.run()
    event = eaux.event()
    lumi  = eaux.luminosityBlock()

    #read all products as specifed in edmCollections
    products = {}
    for k, v in edmCollections.iteritems():
      events.getByLabel(v['label'], v['handle'])
      products[k] = v['handle'].product()

    print run,lumi,event

    for proton in products['ptracks']:

        ismultirp   = -999
        decRPId     = -999
        armId       = -999
        th_y        = -999
        th_x        = -999
        t           = -999
        xi          = -999

        if not proton.valid(): continue

        th_y = (proton.direction().y()) / (proton.direction().mag())
        th_x = (proton.direction().x()) / (proton.direction().mag())
        xi = proton.xi()

        # t
        m = 0.938 # GeV # proton mass
        p = 6500. # GeV # beam energy

        t0 = 2.*m*m + 2.*p*p*(1.-xi) - 2.*sqrt( (m*m + p*p) * (m*m + p*p*(1.-xi)*(1.-xi)) )
        th = sqrt(th_x * th_x + th_y * th_y)
        S = sin(th/2.)
        t = t0 - 4. * p*p * (1.-xi) * S*S

        if proton.method == ROOT.reco.ProtonTrack.rmSingleRP:
            rpId = ROOT.CTPPSDetId( int(proton.contributingRPIds.begin()) )
            decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp()
            ismultirp = 0
        if proton.method == ROOT.reco.ProtonTrack.rmMultiRP:
            rpId = ROOT.CTPPSDetId( int(proton.contributingRPIds.begin()) )
            armId = rpId.arm()
            ismultirp = 1

        reco_protons.append( { 
            'xi' :proton.xi(),
            'thx':th_x,
            'thy':th_y,
            't':t,
            'ismultirp':ismultirp,
            'rpId':rpId
            })
            
    if products['ptracks'].size()>0: break
