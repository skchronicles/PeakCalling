# Fraction of Top N peaks within 100 bp of a STAT1 motif for the 4 methods given

library(ggplot2)

data<-read.table("ALLCALLERS_fractionTOPNPEAKS.txt",header=T)
ggplot(data, aes(x=NPEAK, y=Fraction, fill=Method,  color=Method)) + 
  geom_point(size=0.5) +
  labs(x="Top N Peaks", y="Fraction witin a 100 bp of STAT1 motif")