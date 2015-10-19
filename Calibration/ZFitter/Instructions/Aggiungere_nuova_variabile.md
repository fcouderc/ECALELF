# Stabilire cosa importare
in src/SmearingImporter.cc devi modificare il metodo Import

# Stabilire come la variabile importata deve essere smearata
in src/RooSmearer.cc nel metodo SetSmearedHisto

# Se ti serve di portarti appresso una variabile (es. per smearare (pt1 + pt2) dovevo conoscere pt1 e pt2
vain in interface/ZeeEvent.hh e aggiungila come membro della classe ZeeEvent