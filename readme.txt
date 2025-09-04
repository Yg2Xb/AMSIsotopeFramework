------------------------------READ ME------------------------------
1. selected Be B C N events from 13.5years IHEP New DST /eos/ams/user/z/zuhao/yanzx/ams_data/amsd69n_IHEPQcutDST (by yanzx, exclude all Z = 1 and Z = 2 events, see /eos/ams/user/z/zuhao/yanzx/ams_software/DSTvdev/analysis_amsd69n.C)

2.Pass All events with innerQ > 3.3 && innerQ < 7.7, with RTi,phystrigger,Fiducial Volume、InnerTrk(Hit,Chi2)、basicTOF(build_type < 10, beta>0.4) Cut, Using YiJia Charge, Tuned ToF Beta in AMS SW, RICH beta with all corrections in AMS SW. 

3. specific cuts, tuning and corrections can be found in the Twiki: https://twiki.cern.ch/twiki/bin/view/AMS/Isotopes

4.tree name: saveTree

5. variables for analysis:
    
    double InnerRig(GBL V6)
    double ToFBeta, NaFBeta, AGLBeta; (NaF and AGL are exclusive)
    double ToFEk,   NaFEk,   AGLEk; (Ek means Kinetic Energy per Nuclear, directly converted from beta)
    
    //--------------------------------------------
    ***unsigned int cutStatus***

    Basic Status Bits (0-6):
    Bit 0: Normal Beryllium Nuclei Selection, pass all tracker cuts and bkg cuts, include basicTOF and UpperTOF Q Cut 
    Bit 1: All TOF(geometry, Beta Quality) Cut && TOFBeta > 0
    Bit 2: All RICH(geometry, rich Charge, Beta Quality) Cut && RICH Beta > 0
    Bit 3: Beyond rigidity cutoff 
    Bit 4: Beyond Beta Cutoff using Z/A = 4/7(most strict) with TOF safety factor
    Bit 5: Beyond Beta Cutoff using Z/A = 4/7(most strict) with NaF safety factor
    Bit 6: Beyond Beta Cutoff using Z/A = 4/7(most strict) with AGL safety factor

    Charge Measurement Related Bits (7-26):(all Q cuts are consistent with twiki)
    Bits 7-8: Inner Q cuts details
    - Bit 7: inner Q range cut
    - Bit 8: inner Q RMS cut

    Bit 9: Upper TOF Q cut

    Bits 10-14: L1 Normal cuts details  
    - Bit 10: L1 Q upper limit
    - Bit 11: L1 Q lower limit
    - Bit 12: L1 quality status
    - Bit 13: L1 XY signal
    - Bit 14: L1 Y direction chi-square

    Bits 15-18: L1 Unbiased cuts details
    - Bit 15: Unbiased L1 Q upper limit
    - Bit 16: Unbiased L1 Q lower limit
    - Bit 17: Unbiased L1 quality status  
    - Bit 18: Unbiased L1 XY signal exists

    Bits 19-21: Background cuts details
    - Bit 19: TWIKI background cut
    - Bit 20: No background cut
    - Bit 21: Strict single track cut

    Bits 22-23: TOF cuts details
    - Bit 22: TOF geometry cut
    - Bit 23: TOF beta quality cut

    Bits 24-26: RICH cuts details  
    - Bit 24: RICH geometry cut
    - Bit 25: RICH beta quality cut
    - Bit 26: RICH charge cut

    Usage Examples:
    // Check if specific cut passed
    bool isTrackerCutPassed = (cutStatus & (1 << 0));
    bool isTOFCutPassed = (cutStatus & (1 << 1));
    bool isRICHCutPassed = (cutStatus & (1 << 2));

    // Check combination of cuts
    bool isTrackerAndRICHPassed = (cutStatus & 0x05) == 0x05;

    //--------------------------------------------

    int BkgStudyCut[sel][bkg][det]:
    [sel = 2]: complete Tracker(exclude bkg cut) + ToF or RICH selection, or Loose(remove some beta quality cut) selection
    [bkg = 3]: With Bkg Reduction cuts, or without, or more strict(ntrack==1)
    [det = 3]: ToF, or NaF, or AGL events

    DetectorAndBkgCut and cutStatus can be cross checked

    bool rich_NaF (NaF and AGL are exclusive, rich_NaF ? event pass NaF : pass AGL)

    all other amstreea branches are saved
