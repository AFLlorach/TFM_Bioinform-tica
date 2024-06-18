# TFM
a)	Funcions MKT:
-	“Functions_SFS_aMKT.R”: funcions per realitzar aMKT usant SFS
-	“Functions_Theta_Int_aMKT.R”: funcions per realitzar aMKT a través de l’estima de variabilitat on usem intervals (estimes usant 5 intervals de freqüències del 20%).
-	“Functions_Theta_Stats_aMKT”: funcions per realitzar aMKT a través de l’estima de variabilitat on usem estadístics de variabilitat clàssics (Fujima, Watterson, Fu&Li, fay&Wu i vhf). vhf és un estadístic creat per nosaltres per a altres freqüpencies
Aquestes funcions requereixen de l’script “asymptoticMK_local”  (Haller i Messer, 2017*) i dels paquets d’R “nsl2” i “proto”.
 
  *B.C. Haller, P.W. Messer. (2017). asymptoticMK: A web-based tool for the asymptotic McDonald-Kreitman test. G3: Genes, Genomes, Genetics 7(5), 1569-1575. doi:10.1534/g3.117.039693

  b)	Funcions SLiM:
-	“slim_template”: fitxer amb les instruccions de la simulació (creació de 2 poblacions, divisió en 2 poblacions, possible canvi en el nombre d’individus i final de simulació).
-	“run_construct_slim_conditions_LOCAL”: fitxer amb les diferents condicions (segons si es dóna o no backgrown selection, selecció beneficiosa i canvi de la població efectiva)

c)	Altres funcions:
-	“sum”: script per a sumar les simulacions de cada condició
-	“boostrap_MKT”: script per a realitzar el bootstrap de les dades obtingudes de la suma de les simulacions de cadascuna de les condicions. 
-	“run_plots_MKTa”: script per a obtenir els plots dels resultats

d)	Dades SLiM: conté els fitxers de les dades de les simulacions, de la suma de les simulacions de cada condició i dels Bootstrap obtinguts.

e)	RESULTATS:
-	Taules de resultats “daf.results”, “Theta.int.results” i “Theta.stats.results”.
-	“Plots_MKT_Bootstrap”: arxiu amb tots els plots dels resultats per a cada condició.
