Per fare il fit a un profilo solo

```
#Per fittare tutta l'estensione in y delle Likelihood
root -l -b
.L macro/fitOneProfile.C
fitOneProfile("profile.root","outDir")
```

```
#Per fittare con un range paragonabile all'errore sul fit
./script/fit.sh outProfile.....(senza .root alla fine)
```
