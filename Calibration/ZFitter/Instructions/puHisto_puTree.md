L'idea e' che nel configFile metti gli istogrammi del PU.

A partire da quelli, si creano i root trees del PU e poi vengono usati quelli.

Quindi se hai direttamente i root trees metti quelli nel configFile, altrimenti
parti dai PUHisto e lancia ZFitter con saveRootMacro. 

In /tmp ti trovi i PUTrees: spostali dove ti fa comodo e mettili nel configFile cosi' la volta dopo non li riproduci