// An example of neutrino-nucleon DIS interaction

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // minimal Q2, number of events to generate.
  double Q2min = 25.;
  int nEvent = 100;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event &event = pythia.event;

  // Set up incoming beams, for frame with unequal beam energies.
  // The beams are back-to-back, but with different energies
  pythia.readString("Beams:frameType = 2");
  pythia.readString("Beams:idA = 2212"); // proton
  pythia.settings.parm("Beams:eA", 0.);  // at rest
  pythia.readString("Beams:idB = 12");   // electron neutrino
  pythia.settings.parm("Beams:eB", 1e6); // 1 PeV

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference) and charged current
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  pythia.readString("ProcessLevel:all = on");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("PartonLevel:all = on");
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Initialize.
  pythia.init();

  // Histograms.
  double Wmax = sqrt(4. * 1 * 1e6);
  Hist Qhist("Q [GeV]", 100, 0., 50.);
  Hist Whist("W [GeV]", 100, 0., Wmax);
  Hist xhist("x", 100, 0., 1.);
  Hist yhist("y", 100, 0., 1.);
  Hist pTehist("pT of scattered neutrino [GeV]", 100, 0., 50.);
  Hist pTrhist("pT of radiated parton [GeV]", 100, 0., 50.);
  Hist pTdhist("ratio pT_parton/pT_neutrino", 100, 0., 5.);

  std::ofstream output("output.txt");

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next())
      continue;
    // Four-momenta of proton, neutrino, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 pLeptonIn = event[4].p();
    Vec4 pLeptonOut = event[6].p();
    Vec4 pBosen = pLeptonIn - pLeptonOut;

    // Q2, W2, Bjorken x, y.
    double Q2 = -pBosen.m2Calc();
    double W2 = (pProton + pBosen).m2Calc();
    double x = Q2 / (2. * pProton * pBosen);
    double y = (pProton * pBosen) / (pProton * pLeptonIn);

    // Fill kinematics histograms.
    Qhist.fill(sqrt(Q2));
    Whist.fill(sqrt(W2));
    xhist.fill(x);
    yhist.fill(y);
    pTehist.fill(event[6].pT());

    // pT spectrum of partons being radiated in shower.
    for (int i = 0; i < event.size(); ++i)
      if (event[i].statusAbs() == 43) {
        pTrhist.fill(event[i].pT());
        pTdhist.fill(event[i].pT() / event[6].pT());
      }

    output << "**** Event: " << iEvent << " ****" << endl;
    for (int i = 0; i < event.size(); ++i)
      if (event[i].status() > 0) {
        output << "id: " << event[i].id() << ", energy: " << event[i].e()
               << endl;
      }
  }
  // End of event loop.
  // Statistics and histograms.
  pythia.stat();
  cout << Qhist << Whist << xhist << yhist << pTehist << pTrhist << pTdhist;

  // Done.
  return 0;
}
