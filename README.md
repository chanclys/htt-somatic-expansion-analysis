# HTT Somatic Expansion Analysis - PacBio PureTarget

Comprehensive analysis of Huntington's disease somatic mosaicism using PacBio long-read sequencing with TRGT (Tandem Repeat Genotyper).

## Overview

This project analyzes HTT CAG repeat expansion and genome-wide repeat instability in three samples:
- **Control_11**: Normal healthy control (24-26 CAG, stable)
- **H1**: Huntington's disease patient 1 (inherited 51 CAG, somatic up to 1394 CAG)
- **H3**: Huntington's disease patient 2 (inherited 63 CAG, somatic up to 1078 CAG)

## Key Findings

### HTT Expansion (Main Result)
- **H1**: 53.6x fold expansion (max), somatic burden 29.9 extra CAGs
- **H3**: 46.9x fold expansion (max), somatic burden 43.1 extra CAGs
- Both show extreme somatic mosaicism with reads up to 1000+ CAG repeats

### Genome-Wide Repeat Instability
- **38 repeat loci** analyzed across all samples
- **SCA27B_FGF14 shows massive expansion**: 2.95x-2.60x in both HD patients
- Suggests system-wide repeat instability beyond HTT alone

## Author

Carlos (chanclys) - ccmarinas@ucsd.edu

## License

MIT
