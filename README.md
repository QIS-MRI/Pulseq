# Pulseq
**********
Matlab and Python codes to write a 3D Phase-Cycled bSSFP Cartesian Sequence in the Pulseq Framework.

**********
## Available codes

1) "write_pulseq_PhaseCycledbSSFP.py"

2) "write_pulseq_PhaseCycledbSSFP.m"

## Sequence parameters

| FOV [mm]     | Spatial Resolution | TR [ms] | TE [ms] | RF Excitation Angle [°] | Number of Phase Cycles | BW [Hz] |
|--------------|--------------------|---------|---------|-------------------------|------------------------|---------|
|  384x192x192 |     384x192x64     |  5.02   |   2.51  |       20                |           18           |   868   |

## Diagram
![Diagram of the 3D Phase-cycled bSSFP sequence. RF pulse application followed by 
acquisition with ADC event and balanced gradients. ](seq_diagram.png)

