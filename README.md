# ProbeStation_dataAnalysis
C++ code that creates root files to store the data from the measurements on the prove stations.
It creates two root trees, one to store the raw date and another one that prepares the data to analyze it later.
To run simply write `root main.cpp` on a terminal with `ROOT`.

09/05/25
ficar bins dimensionals pa que el mig del bin "1" caiga damunt de la posicio x=1, utilitzar 0.5 per a ello. Elegir color que tinga molta variacio per els colors mitjans i pa daltonics: https://root.cern.ch/doc/master/classTColor.html#C06a

colores, 0 no es necesario, titulo eje z (colorbar?)

dubte: trac cs del plateau en log scale i pase a escala normal (pF) o faig un fit en la escala normal pa traureu?
colors, automatitzar punts pa fits, valors... preguntar si propagacio de errors be 

automatizar nom cv graphs->guardar en tree nom del file?
IV curves
programacio modular

per a cosaxch millor scatter? tambe serviria per a fer la estimacio de fora de rango del chi2?

iv rango en stv, aplicar a tots, imagen histogramas y canvas 

comparar resultats entre sensors, pasarselos a melissa. 