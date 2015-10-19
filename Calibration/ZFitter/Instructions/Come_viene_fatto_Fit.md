```
./script/Likelihood_fitter.sh
```


Il Fit viene fatto in 

macro/macro_fit.C

Il vero fit viene fatto in IterMinumumFit()

al che Plot()

Gli step in Likelihood_fitter sono:

1) Fare un .root globale con tutti i .root (con haddTGraph)

2) Viene scritto in tmp/ un fitProfile.C

3) Viene lanciato: root -l -b tmp/fitProfile.C


