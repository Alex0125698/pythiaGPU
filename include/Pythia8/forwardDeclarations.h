
namespace Pythia8 {


// I guess Hook -> some external code to do something?

// (Analysis.h)

// performs sphericity analysis on an event
class Sphericity;
// performs thrust analysis on an event
class Thrust;
// helper class to ClusterJet for a jet and its contents
class SingleClusterJet;
// performs a jet clustering
class ClusterJet;
// helper class to CellJet for a cell and its contents
class SingleCell;
// helper class to CellJet for a jet and its contents
class SingleCellJet;
// performs a cone jet search in (eta, phi, E_T) space.
class CellJet;
// Base class, used to derive your own class with your selection criteria
class SlowJetHook;
// helper class to SlowJet for a jet and its contents
class SingleSlowJet;
// performs a recombination jet search in (y, phi, pT) space
class SlowJet;

// (Basics.h)

// the base class for external random number generators
class RndmEngine;
// random number generation according to the Marsaglia-Zaman-Tsang algorithm.
class Rndm;
// four-vectors, in energy-momentum space
class Vec4;
// matrices that encode an arbitrary combination of rotations and boosts
class RotBstMatrix;
// handles a single histogram
class Hist;

// (BeamParticle.h)

// holds info on a parton resolved inside the incoming beam
class ResolvedParton;
// holds info on a beam particle
class BeamParticle;

// (BeamRemnants.h)

// matches the kinematics of the hard-scattering subsystems to that of the two beam remnants
class BeamRemnants;

// (BeamShape.h)

// Base class to set beam momentum and interaction spot spread
class BeamShape;

// (BoseEinstein.h)

// a simple container for studied hadrons
class BoseEinsteinHadron;
// shifts the momenta of identical particles relative to each other, to simulate Bose-Einstein effects to some approximation
class BoseEinstein;

// (ColourReconnection.h)

// Contain a single colour chain. It always start from a quark and goes to an anti quark or from an anti-junction to a junction
class ColourDipole;
// Junction class. In addition to the normal junction class, also contains a list of dipoles connected to it.
class ColourJunction; // : public Junction
// ??
class TrialReconnection;
// ??
class ColourParticle; // : public Particle
// handles the colour reconnection
class ColourReconnection;

// (colourTracing.h)

// used to trace colours within the event record
class ColourTracing;

// (Event.h)

// holds info on a particle in general
class Particle; 
// stores what kind of junction it is, the colour indices of the legs at the junction and as far out as legs have been traced
class Junction;
// holds all info on the generated event
class Event;

// (FJcore.h)

// skip

// (FragmentationFlavZpT.h)

// simple container for flavour, including the extra properties needed for popcorn baryon handling
class FlavContainer;
// used to select quark and hadron flavours
class StringFlav;
// used to sample the fragmentation function f(z)
class StringZ;
// used to select select transverse momenta
class StringPT;

// (FragmentationSystems.h)

// ColSinglet class contains info on an individual singlet. Only to be used inside ColConfig, so no private members.
class ColSinglet;
// describes the colour configuration of the whole event
class ColConfig;
// contains the information related to one string section in the evolution of a multiparton system. Only to be used inside StringFragmentation and MiniStringFragmentation
class StringRegion;
// contains the complete set of all string regions. Only to be used inside StringFragmentation, so no private members
class StringSystem;

// (HadronLevel.h)

// contains the top-level routines to generate the transition from the partonic to the hadronic stage of an event
class HadronLevel;

// (HadronScatter.h)

// ??
class SigmaPartialWave;
// hold details of a pair of hadrons which will scatter. Stores indices in event record and the measure used for ordering
class HadronScatterPair;
// ??
class HadronScatter;

// (HardDiffraction.h)

// handles hard diffraction, together with PartonLevel
class HardDiffraction;

// (HelicityBasics.h)

// complex four-vector wave functions. The Wave4 class can be multiplied with the GammaMatrix class to allow for the writing of helicity matrix elements
class Wave4;
// a special sparse matrix class used to write helicity matrix elements in conjuction with the Wave4 class
class GammaMatrix;
// containing helicity information
class HelicityParticle; // : public Particle

// (HelicityMatrixElements.h)

// used in tau tacays
class HelicityMatrixElement;

// for the hard process of two fermions -> W/W' -> two fermions.
class HMETwoFermions2W2TwoFermions; // : public HelicityMatrixElement;
// for the hard process of two fermions ->  photon/Z/Z' -> two fermions.
class HMETwoFermions2GammaZ2TwoFermions; // : public HelicityMatrixElement;
// for the decay of a Higgs ->  two fermions
class HMEHiggs2TwoFermions; // : public HelicityMatrixElement;
// for all tau decay helicity matrix elements
class HMETauDecay; // : public HelicityMatrixElement;
// for the hard process of X -> two fermions
class HMEX2TwoFermions; // : public HelicityMatrixElement;

// for the hard process of W/W' -> two fermions
class HMEW2TwoFermions; // : public HMEX2TwoFermions;
// for the hard process of photon -> two fermions
class HMEGamma2TwoFermions; // : public HMEX2TwoFermions;
// for the hard process of Z/Z' -> two fermions
class HMEZ2TwoFermions; // : public HMEX2TwoFermions;

class HMETau2Meson; // : public HMETauDecay;
class HMETau2TwoLeptons; // : public HMETauDecay;
class HMETau2TwoMesonsViaVector; // : public HMETauDecay;
class HMETau2TwoMesonsViaVectorScalar; // : public HMETauDecay;
class HMETau2ThreeMesons; // : public HMETauDecay;
class HMETau2ThreePions; // : public HMETau2ThreeMesons;
class HMETau2ThreeMesonsWithKaons; // : public HMETau2ThreeMesons;
class HMETau2ThreeMesonsGeneric; // : public HMETau2ThreeMesons;
class HMETau2TwoPionsGamma; // : public HMETauDecay;
class HMETau2FourPions; // : public HMETauDecay;
class HMETau2FivePions; // : public HMETauDecay;
class HMETau2PhaseSpace; // : public HMETauDecay;

// (HiddenValleyFragmentation.h)

// used to select HV-quark and HV-hadron flavours.
class HVStringFlav; // : public StringFlav;
// used to select select HV transverse momenta.
class HVStringPT; // : public StringPT;
// used to sample the HV fragmentation function f(z).
class HVStringZ; // : public StringZ;
// routines to fragment a Hidden Valley partonic system
class HiddenValleyFragmentation;

// (History.h)

// holds information about one radiator, recoiler, emitted system. This class is a container class for History class use
class Clustering;
// A History object represents an event in a given step in the CKKW-L
// clustering procedure. It defines a tree-like recursive structure,
// where the root node represents the state with n jets as given by
// the matrix element generator, and is characterized by the member
// variable mother being null. The leaves on the tree corresponds to a
// fully clustered paths where the original n-jets has been clustered
// down to the Born-level state. Also states which cannot be clustered
// down to the Born-level are possible - these will be called
// incomplete. The leaves are characterized by the vector of children
// being empty.
class History;

// (Info.h)

// The Info class contains a mixed bag of information on the event
// generation activity, especially on the current subprocess properties,
// and on the number of errors encountered.
class Info;

// (JunctionSplitting.h)

// JunctionSplitting takes an event and seperate junction chains from each other, such that no junctions are colour connected to each other
class JunctionSplitting;

// (LesHouches.h)

// processes stored in LHAup
class LHAProcess;
// particles stored in LHAup
class LHAParticle;
// base class for initialization and event information from an external parton-level generator
class LHAup;
// information read from a Les Houches Event File
class LHAupLHEF; // : public LHAup;
// information read from PYTHIA 8 itself, for output
class LHAupFromPYTHIA8; // : public LHAup;
// with LHEF 3.0 information read from PYTHIA 8 itself, for output.
class LHEF3FromPythia8; // : public LHAup;

// (LHEF3.h) Les Houches Event File version 3

// The XMLTag struct is used to represent all information within an XML tag.
// It contains the attributes as a map, any sub-tags as a vector of pointers
// to other XMLTag objects, and any other information as a single string.
struct XMLTag;
// represents the information in a weights tag
struct LHAweights;
// Collect different scales relevant for an event
struct LHAscales;
// Collect generator information for an event file
struct LHAgenerator;
// Collect the wgt information
struct LHAwgt;
// Collect the weight information
struct LHAweight;
// assigns a group-name to a set of LHAweight objects
struct LHAweightgroup;
// assigns a group-name to a set of LHAwgt objects
struct LHArwgt;
// assigns a group-name to a set of LHAweightgroup objects
struct LHAinitrwgt;
// The HEPRUP class is a simple container corresponding to the Les Houches
// accord common block with the same name. The members are named in the same
// way as in the common block.
class HEPRUP;
// The HEPEUP class is a simple container corresponding to the Les Houches
// accord common block with the same name. The members are named in the same
// way as in the common block.
class HEPEUP;
// The Reader class is initialized with a stream from which to read a
// version 1/3 Les Houches Accord event file. In the constructor of
// the Reader object the optional header information is read and then
// the mandatory init is read. After this the whole header block
// including the enclosing lines with tags are available in the public
// headerBlock member variable. Also the information from the init
// block is available in the heprup member variable and any additional
// comment lines are available in initComments. After each successful
// call to the readEvent() function the standard Les Houches Accord
// information about the event is available in the hepeup member
// variable and any additional comments in the eventComments
// variable. A typical reading sequence would look as follows:
class Reader;
// The Writer class is initialized with a stream to which to write a
// version 1.0 or 3.0 Les Houches Accord event file. In the init() function of
// the Writer object the main XML tag, header and init blocks are written,
// with the corresponding end tag is written by print_end_tag().
// After a Writer object (in the following called "writer") has been created,
// it is possible to assign version (3 by default) information.
class Writer;

// (Merging.h)

// wrapper class for the interface of matrix element merging and Pythia8.
class Merging;

// (MergingHooks.h)

// This class holds information on the desired hard 2->2 process
// for the merging. This class is a container class for History class use.
class HardProcess;
// base class for user input to the merging procedure
class MergingHooks;

// (MiniStringFragmentation.h)

// The MiniStringFragmentation class contains the routines to fragment
// occasional low-mass colour singlet partonic systems, where the string
// approach is not directly applicable (for technical reasons).
class MiniStringFragmentation;

// (MultipartonInteractions.h)

// SigmaMultiparton is a helper class to MultipartonInteractions.
// It packs pointers to the allowed processes for different
// flavour combinations and levels of ambition.
class SigmaMultiparton;
// The MultipartonInteractions class contains the main methods for the
// generation of multiparton parton-parton interactions in hadronic collisions.
class MultipartonInteractions;

// (ParticleData.h)

// holds info on a single decay channel
class DecayChannel;
// holds info on a single particle species
class ParticleDataEntry;
// holds a map of all ParticleDataEntries
class ParticleData;

// (ParticleDecays.h)

// base class for the external handling of decays
class DecayHandler;
// contains the routines to decay a particle
class ParticleDecays;

// (PartonDistributions.h)

// Base class for parton distribution functions
class PDF;
// different types of PDFs
class GRV94L; // : public PDF;
class CTEQ5L; // : public PDF;
class MSTWpdf; // : public PDF;
class CTEQ6pdf; // : public PDF;
class ProtonPoint; // : public PDF;
class GRVpiL; // : public PDF;
class PomFix; // : public PDF;
class PomH1FitAB; // : public PDF;
class PomH1Jets; // : public PDF;
class NNPDF; // : public PDF;
class LHAPDF; // : public PDF;
// Gives electron (or muon, or tau) parton distribution
class Lepton; // : public PDF;
// Gives electron (or other lepton) parton distribution when unresolved
class LeptonPoint; // : public PDF;
// Gives neutrino parton distribution when unresolved (only choice for now)
class NeutrinoPoint; // : public PDF;

// (PartonLevel.h)

// contains the top-level routines to generate the partonic activity of an event
class PartonLevel;

// (PartonSystems.h)

// contains info on an individual singlet
class PartonSystem;
// describes the whole set of subcollisions
class PartonSystems;

// (PhaseSpace.h)

// base class for  phase space generators used in the selection of hard-process kinematics
class PhaseSpace;
// derived classes for each type of final state
class PhaseSpace2to1tauy; // : public PhaseSpace;
class PhaseSpace2to2tauyz; // : public PhaseSpace;
class PhaseSpace2to2elastic; // : public PhaseSpace;
class PhaseSpace2to2diffractive; // : public PhaseSpace;
class PhaseSpace2to3diffractive; // : public PhaseSpace;
class PhaseSpace2to2nondiffractive; // : public PhaseSpace;
class PhaseSpace2to3tauycyl; // : public PhaseSpace;
class PhaseSpace2to3yyycyl; // : public PhaseSpace;
// Les Houches events
class PhaseSpaceLHA; // : public PhaseSpace;

// (ProcessContainer.h)

// The ProcessContainer class combines pointers to matrix element and
// phase space generator with general generation info.
class ProcessContainer;
// The SetupContainers class turns the list of user-requested processes
// into a vector of ProcessContainer objects, each with a process.
class SetupContainers;

// (ProcessLevel.h)

// The ProcessLevel class contains the top-level routines to generate
// the characteristic "hard" process of an event.
class ProcessLevel;

// (Pythia.h)

// contains the top-level routines to generate an event.
class Pythia;

// (PythiaComplex.h)

// none

// (PythiaStdlib.h)

// none

// (ResonanceDecays.h)

// The ResonanceDecays class handles the sequential decay of resonances
// that are part of the hard process (t, W, Z, H, SUSY,...).
class ResonanceDecays;

// (ResonanceWidths.h)

// the base class. Also used for generic resonaces
class ResonanceGeneric; // : public ResonanceWidths;
// derived classes for specific resonances
class ResonanceGmZ; // : public ResonanceWidths;
class ResonanceW; // : public ResonanceWidths;
class ResonanceTop; // : public ResonanceWidths;
class ResonanceFour; // : public ResonanceWidths;
class ResonanceH; // : public ResonanceWidths;
class ResonanceHchg; // : public ResonanceWidths;
class ResonanceZprime; // : public ResonanceWidths;
class ResonanceWprime; // : public ResonanceWidths;
class ResonanceRhorizontal; // : public ResonanceWidths;
class ResonanceExcited; // : public ResonanceWidths;
class ResonanceGraviton; // : public ResonanceWidths;
class ResonanceKKgluon; // : public ResonanceWidths;
class ResonanceLeptoquark; // : public ResonanceWidths;
class ResonanceNuRight; // : public ResonanceWidths;
class ResonanceZRight; // : public ResonanceWidths;
class ResonanceWRight; // : public ResonanceWidths;
class ResonanceHchgchgLeft; // : public ResonanceWidths;
class ResonanceHchgchgRight; // : public ResonanceWidths;

// (RHadrons.h)

// The RHadrons class contains the routines for the production and decay
// of long-lived heavy coloured particles.
class RHadrons;

// (Settings.h)

// classes for each type of setting
class Flag;
class Mode;
class Parm;
class Word;

// classes for arrays of settings
class FVec;
class MVec;
class PVec;

// This class holds info on flags (bool), modes (int), parms (double),
// words (string), fvecs (vector of bool), mvecs (vector of int) and pvecs
// (vector of double).
class Settings;

// (SigmaProcess.h)

// helper class for partons and their flux in a beam
class InBeam;
// helper class for colliding parton pairs and their flux
class InPair;
// base class for cross section calculations
class SigmaProcess;
// base class for unresolved and minimum-bias processes
class Sigma0Process; // : public SigmaProcess;
// base class for 2 -> 1 processes
class Sigma1Process; // : public SigmaProcess;
// base class for 2 -> 2 processes
class Sigma2Process; // : public SigmaProcess;
// base class for 2 -> 3 processes
class Sigma3Process; // : public SigmaProcess;
// wrapper class for Les Houches Accord external input
class SigmaLHAProcess; // : public SigmaProcess;

// (SigmaTotal.h)

// The SigmaTotal class contains parametrizations of total, elastic and
// diffractive cross sections, and of the respective slope parameter.
class SigmaTotal;

// (SigmaCompositeness.h)
// (SigmaEW.h)
// (SigmaExtraDim.h)
// (SigmaGeneric.h)
// (SigmaHiggs.h)
// (SigmaLeftRightSym.h)
// (SigmaLeptoquark.h)
// (SigmaNewGaugeBosons.h)
// (SigmaOnia.h)
// (SigmaQCD.h)
// (SigmaSUSY.h)

// (SLHAinterface.h)

// The SLHAinterface class handles communication between Pythia and SusyLesHouches.
class SLHAinterface;

// (SpaceShower.h)

// Data on radiating dipole ends, only used inside SpaceShower
class SpaceDipoleEnd;
// does spacelike showers
class SpaceShower;

// (StandardModel.h)

// The AlphaStrong class calculates the alpha_strong value at an arbitrary
// scale, given the value at m_Z, to zeroth, first or second order.
class AlphaStrong;
// The AlphaEM class calculates the alpha_electromagnetic value at an
// arbitrary scale, given the value at 0 and m_Z, to zeroth or first order.
class AlphaEM;
// The CoupSM class stores and returns electroweak couplings,
// including Cabibbo-Kobayashi-Maskawa mass mixing matrix elements.
class CoupSM;
// Generic couplings class base class
class Couplings; // : public CoupSM;

// (Streams.h)

// skip

// (StringFragmentation.h)

// The StringEnd class contains the information related to
// one of the current endpoints of the string system.
// Only to be used inside StringFragmentation, so no private members.
class StringEnd;
// The StringFragmentation class contains the top-level routines
// to fragment a colour singlet partonic system.
class StringFragmentation;

// (StringLength.h)

// used to calculate the lambda measure of strings and junctions
class StringLength;

// (SusyCouplings.h)

// Auxiliary class to compute and store various SM and SUSY couplings.
class CoupSUSY; // : public Couplings;

// (SusyLesHouches.h)

// the generic SLHA block
template <class T> class LHblock;
// Derived class for generic blocks containing vectors of strings
class LHgenericBlock; // : public LHblock<string>;
// the generic SLHA matrix block
template <int size> class LHmatrixBlock;
// generic SLHA tensor block
template <int size> class LHtensor3Block;
// ??
class LHdecayChannel;
// ??
class LHdecayTable;
// ??
class SusyLesHouches;

// (SusyResonanceWidths.h)

// more resonance widths for susy particles
class SUSYResonanceWidths; // : public ResonanceWidths;
class ResonanceSquark; // : public SUSYResonanceWidths;
class ResonanceGluino; // : public SUSYResonanceWidths;
class ResonanceNeut; // : public SUSYResonanceWidths;
class ResonanceChar; // : public SUSYResonanceWidths;
class ResonanceSlepton; // : public SUSYResonanceWidths;

// (SusyWidthFunctions.h)

// base class for SUSY 3-body decay width functions.
class WidthFunction;
// ??
class StauWidths; // : public WidthFunction;

// (TauDecays.h)

// This class decays tau leptons, with helicity information
class TauDecays;

// (TimeShower.h)

// Data on radiating dipole ends; only used inside TimeShower class.
class TimeDipoleEnd;
// does timelike showers
class TimeShower;

// (UserHooks.h)

// to allow user access to program at different stages
class UserHooks;
// SuppressSmallPT is a derived class for user access to program execution.
// It is a simple example, illustrating how to suppress the cross section
// of 2 -> 2 processes by a factor pT^4 / (pT0^2 + pT^2)^2, with pT0 input,
// and also modify alpha_strong scale similarly.
class SuppressSmallPT; // : public UserHooks;

// (WeakShowerMEs.h)

// provides ME's needed for W/Z emission in ISR or FSR
class WeakShowerMEs;

}