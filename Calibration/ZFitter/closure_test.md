1) Avere un initFile con scala e sigma da iniettare ai Toy

2) Solo i file con label "s" vengono utilizzati nel closure test: quindi solo il MC

3) In SmearingImporter si decide la divisione in dati vs MC (in base a %5)

3) in src/RooSmearer.cc Init puoi decidere come inizializzare la scala/sigma di quello che hai chiamato MC

```
if(mcToy && false){

do something; puoi randomizzare ad esempio la scala e la sigma ecc...
}
```