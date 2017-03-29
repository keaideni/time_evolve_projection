#!/bin/sh

#mkdir Corr
mkdir data
mkdir result

cat > data/QNosave.txt <<EOF
LatticeSize= 4
ParticleNo= 2
SiteNo= 2
DeltaQL= 2
DeltaQR= 3
D= 200
SweepNo= 3
EdgeCondition= 1
Energy= 0
Wz= 1
Wc= 1
gr= 0.0055
gl= 0.0245
EOF
