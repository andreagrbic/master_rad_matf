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
a) Aritmetička sredina cijena akcija
b) Evropska kol opcija čija je cijena eksplicitno data Blek-Šolsovom formulom
c) Geometrijska azijska opcija za čiju cijenu postoji eksplicitna formula
