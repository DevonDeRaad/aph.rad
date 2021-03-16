# MSG Illumina Deep Sequencing Mapping Protocol

1. Isolate genomic DNA using SPRI bead protocol: <https://github.com/DevonDeRaad/aph.rad/blob/master/dna.extraction.protocol.md>
2. DNA Quantification using Qubit dsDNA Assay - 1&micro;L
3. Dilute DNA from each sample to a standard concentration DNA – dilute with Epure H20 – (10ng/&micro;L)

### Digest the genomic DNA samples with NdeI.
Prepare a digestion master mix. Add per sample:	
- 2 &micro;L NEB CutSmart Buffer
- 0.15 &micro;L NdeI (3U) use NEB enzyme conc : 4000U of 20,000U/&micro;L 
- 7.85 &micro;L water
---
- 10 &micro;L total (per sample)

For each sample: add 10 &micro;L digestion master mix to 10 &micro;L sample (100ng DNA)
- Incubate at 37 &deg;C for 3h in a thermocycler. 
- Inactivate at 65 &deg;C for 20 minutes in a thermocycler.

### Ligate bar-coded adapters and pool samples.
Make ligation mastermix. Add per sample:
- 5 µL10X T4 DNA Ligase Buffer
- 0.38 &micro;L T4 ligase (150 U) use NEB enzyme conc : 100,000U of 400,000U/mL
- 23.62 &micro;L water
---
-29 &micro;L total (per sample)
 
- Use plate of SJM 40 uM annealed adapters that were diluted to 5 uM 
- Add to each digested sample 1 &micro;L of adapter oligos (5 uM stock solns), then 29 &micro;L of ligation mastermix:
- Ligate in a thermocycler at 25 &deg;C for 2 hours. Inactivate at 65 &deg;C for 10 minutes.

### Isopropanol precipitation
- Use pipette to take up 50 &micro;L from each well – combine all liquid into a 15 mL tube
- Add: .1x sample volume of 3M sodium acetate pH5.2 - e.g. 100 samples > 500 &micro;L
- Add: 1x sample volume of isopropanol - 50 &micro;L - e.g. 100 samples > 5000 &micro;L
- Mix the 15mL tube well by inverting several times
- Aliquot the sample evenly between 4; 1.7mL eppi tubes
- Chill overnight at 4 &deg;C for isopropanol.
- Pellet the precipitate by centrifuging at 14000X for 9 min.
- Carefully pour off the supernatant.  
- Wash the sides of the tube with cold 70% ethanol.  
- Centrifuge at 4000X for another 5 minutes and remove the supernatant.  
- Speedvac for 1-2 min on medium heat
- Resuspend the pellet in each 1.5 mL tube in 25 &micro;L 10mM Tris-HCl + 0.1% Tween20, pH 8.5.  
- Heat 30 minutes at 65 &deg;C.  
- Pool the 4 1.5 mL Eppis from each plate in new tube – end with ~100 &micro;L per sample

### Phenol:Chloroform Extraction
- Add 1x sample volume of phenol:chloroform - 100 &micro;L
- Mix gently; spin for 10 minutes at max speed.
- Remove aqueous layer (top layer) to a new tube. Dump the bottom layer into phenol:chloroform liquid waste container.
- Add 1x remaining sample volume (after removing aqueous layer) of chloroform.
- Mix; spin 5 minutes at max speed.
- Remove aqueous (top) layer to a new tube. Dump bottom layer into waste container.

record volumes:
| Sample name | Sample volume | P:C volume | Sample volume | Chloroform volume | date completed |
|-------------|---------------|------------|---------------|-------------------|----------------|
|-            |-              |-           |-              |-                  |-               |

### Bead purify using the Agencourt AMPure PCR purification kit.  	
- Swirl bottle to resuspend beads.
- Check volume of the sample
- Add 1.5x sample volume of beads, mix.
- Incubate at room temp for 5 min.
- Place on magnet for 10 min.
- Aspirate and discard supernatant.
- Add 200 &micro;L of fresh ETOH 80%.
- Incubate at room temp for 1 min.
- Aspirate and discard supernatant.
- Repeat ETOH wash.
- Take off magnet - Dry at 37 &deg;C for 2-3 minutes
- Off the magnet, resuspend in 32 &micro;L of 10mM Tris-HCl + 0.1% Tween20, pH 8.5.
- Incubate at room temp 1 min.
- Place back on magnet for 5 min.
- Transfer supernatant to a new tube - Spec 1 &micro;L of DNA on Qubit.

record sample volume post Phenol:Chloroform & Qubit results post bead purification
| Sample name | Sample volume | bead volume | Qbit conc. ng/uL | uL DNA | total ng |
|-------------|---------------|-------------|------------------|--------|----------|
|-            |-              |-            |-                 |-       |-         |

### Run pooled sample on a BluePippin Prep 2% Cassette with internal standards.
- Use cassette 2% DF Marker V1
- Make sure pooled sample is at 30&micro;L – otherwise bring volume up to 30&micro;L with 1XTE
- Allow loading solution and sample to equilibrate to RT
- Combine 30&micro;L DNA sample with 10&micro;L loading solution/marker mix (labeled marker V1)
- Run an internal marker lane
- Mix samples and briefly centrifuge
- Set Pippin Prep to elute the 495-605 bp range of the DNA fragments

### Buffer Exchange sample following elution from Pippin cassette
- Swirl bottle to resuspend beads.
- Check volume of the sample
- Add 2X sample volume of beads, mix.
- Incubate at room temp 5 min.
- Place on magnet 10 min.
- Aspirate and discard supernatant.
- Add 200 µl of fresh ETOH 80%.
- Incubate at room temp 1 min.
- Aspirate and discard supernatant.
- Repeat ETOH wash.
- Take off magnet - Dry at 37 C for 3 minutes
- Off the magnet, resuspend in 20 µl of 10mM Tris-HCl + 0.1% Tween20, pH 8.5
- Incubate at room temp 1 min.
- Place back on magnet for 5 min.
- Transfer supernatant to a new tube.
- Spec 1 &micro;L of DNA on Qubit. 

Record volumes post pippin prep & Qubit post pippin + buffer exchange
| Sample name | Sample volume | bead volume | Qbit conc. ng/uL | uL DNA | total ng |
|-------------|---------------|-------------|------------------|--------|----------|
|-            |-              |-            |-                 |-       |-         |

### Amplify bar-coded fragments using the Phusion PCR kit and FC1 and FC2 primers 
Make 2 50 &micro;L PCR reactions for each plate with the following:
- 10 &micro;L 5X Phusion buffer
- 1.0 &micro;L 10mM dNTP
- 2.5 &micro;L FC1 (10 µM) primer
- 2.5 &micro;L FC2 (10 µM) primer
- 0.5 &micro;L Phusion enzyme
- x &micro;L (x = 1-2 ng template DNA) ~ use 2ng of template if possible 
- 33.5 - x &micro;L water

Run the FC_PCR program for 13 cycles
The program should be:
1. 98 &deg;C 30 sec
2. 98 &deg;C 10 sec ----
3. 62 &deg;C 15 sec     13 cycles
4. 72 &deg;C 15 sec ----
5. 72 &deg;C 7 min
6. 4 &deg;C	        -	
- Move PCR product to 1.7mL tube
- Save 1 &micro;L of PCR product to Qubit

Record PCR setup & needed bead purification amounts
| Sample name | conc. ng/uL | uL DNA | uL Epure H20 | FC1 primer | volume post PCR | bead volume (.8x) |
|-------------|-------------|--------|--------------|------------|-----------------|-------------------|
|-            |-            |-       |-             |-           |-                |-                  |

### Bead purify using AMPure Beads to remove all primer dimers and DNA up to 150 bp  
- Combine all PCR reactions from a single plate into 1 bead clean up tube
- Check volume of each sample
- Swirl bottle to resuspend beads.
- Add (0.8x sample volume) of beads, mix.
- Incubate at room temp 5 min.
- Place on magnet 10 min.
- Aspirate and discard supernatant – Save 1 &micro;L sup to Qubit
- Add 200 &micro;L of fresh ETOH 80%.
- Incubate at room temp 1 min.
- Aspirate and discard supernatant.
- Repeat ETOH wash.
- Take tubes off magnet - Dry at 37 C for 3 mins
- Off the magnet, resuspend in 15 &micro;L of 10mM Tris-HCl + 0.1% Tween20, pH 8.5
- Incubate at room temp 1 min.
- Place back on magnet for 5 min.
- Transfer supernatant to a new tube.
- Spec 1 &micro;L of Final MSG Library on Qubit.

| Sample name | post PCR conc. ng/uL | bead discard conc. ng/uL | final MSG library conc. ng/uL | uL left  | bp size  | total nanograms | FC1 primer |
|-------------|----------------------|--------------------------|-------------------------------|----------|----------|-----------------|------------|
|-            |-                     |-                         |-                              |-         |-         |-                |-           |