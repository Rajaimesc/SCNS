From the CAFANA github repo (https://github.com/cafana) :
"CAFAna is a set of C++ tools designed to facilitate construction and analysis of distributions from event records in high energy physics datasets. Its main design goals are speed, modularity, and simplicity and consistency of user experience. It evolved as the main day-to-day analysis framework used by the NOvA experiment. CAFAna provides tools for:

fast iteration over data files to produce histograms
extensive efficient histogram manipulation (using Eigen as a backend)
a library with a number of "calculators" for neutrino oscillation probabilities in various contexts
means for efficient fitting of histograms produced by the framework under neutrino oscillation or systematic uncertainty variations.
The name originates from "Common Analysis Files" (or "Common Analysis Format" files; the oral history is conflicted)—CAFs—which remain the foundation of the approach, but the system has since been generalized beyond the NOvA format. Current applications include efforts on NOvA as well as users from the DUNE collaboration and the Fermilab SBN program."

The application of CAFANA to the short baseline neutrino program is done with the addition of the library SBNAna, which contains the declarations of the variables reconstrcuted from the DAQ in each detector (https://github.com/SBNSoftware/sbnanaobj/tree/89211c82e163f60752afbae6526787af3a317437/sbnanaobj/StandardRecord), and the Core(https://github.com/cafana/CAFAnaCore) designed to handle data in an histogram format.

ICARUS consists of a time projection chamber filled with Liquid Argon designed to collect inverse beta decays products from the muon neutrino beams of NuMI and BNB located at the fermilab. This with the aim of measuring eventual electronic beta decays in the detector that will indicate an oscillation process that for the distnances at which the detector is placed will indicate evidence that backs the existence of the sterile neutrino.
Being located at ground level ICARUS receives a huge amount of cosmic rays for every beam run in whith it does DAQ, for this region a subdetector specialized for the cosmic rays has been built on top of it. The three subdetectors of ICARUS are then:

1. A time projection chamber that collect the ionized particles created by the drift of an incoming particle from the beam or cosmic background that reconstructs the position and topology of the event.
2. A photomultiplier system that collects the time information of the event along the light produced in this
3. A cosmic ray tagger outside the fiducial volume that triggers a signal whenever a cosmic ray crosses it, allowing an association with the events that occur within the TPC.

![image](https://github.com/Rajaimesc/SCNS/assets/20934233/88d41258-ac2b-4439-a92c-5d09aa4d790e) #A schematic view of ICARUS

The reconstruction of the tracks is done via an algorithm called Pandora that classifies and discriminates tracks as cosmic rays or possible neutrino events. However, the efficiency of Pandora in this task is not the best and a further classification needs to be done to discriminate possible signal from cosmic rays with the recontructed variables used. This process is called TripleMatch and consists in the association of the three signals mentioned before. The algorithm will be shown stepby step in the two files TripleMatch.h and TripleMatch.C
Each call and function explained step by step along the syntaxis of CAFE (CAFana Compilator)



