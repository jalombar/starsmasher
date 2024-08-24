#!/bin/sh

# Clean
rm DENS.txt TIME.txt DENS_TIME.txt

# This is the time, converted to seconds, conversion factor = sqrt(Rsun^3/(G*Msun)))
cat energy0.sph| awk '{print $1*1593.0}' >> TIME.txt

# This is rho, converted to g/cm^3, conversion factor = Msun/Rsun^3
cat log0.sph  | grep "rhomax= " | awk '{print $2*5.91}'  >> DENS.txt

# Paste them
paste TIME.txt DENS.txt | column -s $'\t' -t > DENS_TIME.txt

# Add a header
sed -i '1i #1: Time (secs) #2: Density (g/cm^3)' DENS_TIME.txt

echo "Remove lines which only have one column in DENS_TIME.txt"
