# Vrednovanje azijskih opcija
Ovdje se nalaze svi iskomentarisani kodovi implemetirani za potrebe master rada i primjeri. 
Za pokretanje kodova od dodatnih biblioteka je potrebna samo biblioteka *MASS* i njena funkcija *rlm* za primjenu robusne linearne regresije. Sve ostalo je ručno implementirano.

# Opis fajlova i funkcija u njima
1. U fajlu GBM.R se nalazi vektorizovana funkcija koja služi za generisanje trajektorija geometrijskog Braunovog kretanja.
2. U fajlu naivni_monte_karlo.R je implementirana ,,naivna'' Monte-Karlo metoda za vrednovanje azijskih opcija. 
3. U fajlu AV_GBM.R se nalazi vektorizovana funkcija kojom se generišu trajektorije geometrijskog Braunovog kretanja koristeći antitetičko uzorkovanje. 
Na kraju je dat primjer generisanja putanja ovom metodom.
4. U fajlu AV_GBM_MC.R se nalazi funkcija za određivanje Monte-Karlo ocjene cijene aritmetičke azijske akcije koja koristi metod antitetičkog uzorkovanja za redukciju disperzije. 
5. U fajlu CV_GBM_MC.R je implementirana metoda kontrolnih promjenljivih za redukciju disperzije Monte-Karlo ocjene. Funkcija je prilagođena trima kontrolnim promjenljivim, i to:
  * Aritmetička sredina cijena akcija
  * Evropska kol opcija čija je cijena eksplicitno data Blek-Šolsovom formulom
  * Geometrijska azijska opcija za čiju cijenu postoji eksplicitna formula
6. U fajlu CV_GBM_MC.R se nalazi funkcija kojom je implementirana kombinovana metoda kontrolnih promjenljivih. Koriste se tri prethodno navedene kontrolne promjenljive.
7. Fajl kov_poredjenje.R služi za poređenje kovarijacija niza ocjena dobijenog ,,naivnom'' Monte-Karlo metodom i niza dobijenog metodom kontrolnih promjenljivih. Uloga posmatranja te kovarijacije jeste određivanje najefikasnije kontrolne promjenljive za redukciju disperzije i smanjenje greške Monte-Karlo ocjene.
8. Fajl vorst_aproksimacija.R sadrži funkcije za računanje Vorstove aproksimacije za cijenu aritmetičke azijske opcije, kao i funkcije za određivanje gornje i donje granice te cijene. 
