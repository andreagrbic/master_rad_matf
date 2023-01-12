# Metoda kontrolnih slučajnih veličina (Control Variates) za redukciju 
# disperzije
# Koristićemo tri tipa kontrolnih slučajnih veličina:
#1) Aritmetičkа sredina cijena akcija.
#2) Evropska kol opcija čija je cijena eksplicitno data B-Š formulom.
#3) Geometrijska azijska opcija za čiju cijenu postoji eksplicitna formula.
# Očekivanja svih kontrolnih promjenljivih pri mjeri neutralnoj od rizika
# su poznate vrijednosti, a ove kontrolne slučajne promjenljive su korelisane
# sa nepoznatom cijenom aritmetičke azijske opcije čiju vrijednost 
# želimo da ocijenimo. 

# Potrebno je navesti početne vrijednosti parametara na početku, jer je prije 
# pokretanja glavne funkcije AzijskaOpcija_CV potrebno izračunati matematička
# očekivanja kontrolnih slučajnih veličina za zadate parametre.

# S_0 - početna cijena akcije
# r - bezrizična kamatna stopa
# sigma - volatilnost akcije
# n - broj vremenskih intervala kojima diskretizujemo slučajni proces, odnosno
# broj koraka
# N - broj trajektorija geometrijskog Braunovog kretanja 
# delta_t - dužina vremenskih intervala kojima diskretizujemo proces neutralan
# od rizika, vremenski intervali su porazumijevano ekvidistantni
# beta - nivo povjerenja za konstrukciju intervala povjerenja
# br_kontrolnih_simulacija - broj "pilot" simulacija kojima nalazimo optimalni
# kontrolni parametar c


library(MASS) 
S_0 <- 100
K <- 70
r <- 0.05
sigma <- 0.3
T <- 1
n <- 100
N <- 10000 
delta_t <- T / n
beta <- 0.95
br_kontrolnih_simulacija <- 1000

#' @title Funkcija za određivanje prve kontrolne slučajne veličine.
#' 
#' @description Računa se prosjek cijena akcija koje prate geometrijsko
#' Braunovo kretanje.
#'
#' @param putanja numeric. Vektor koji predstavlja jednu putanju geometrijskog
#' Braunovog kretanja - jednu putanju promjena kretanja cijena akcija.
#'
#' @return prosjek cijena akcija za datu putanju/simulaciju
KontrolnaTip1 <- function(putanja) {
  mean(putanja[-1])
}

# E_KontrolnaTip1 - matematičko očekivanje prve kontrolne slučajne veličine u 
# odnosu na mjeru neutralnu od rizika.

E_KontrolnaTip1 <-
  (S_0 / n) * exp(r * delta_t) * (1 - exp(r * delta_t * n)) / (1 - exp(r * delta_t))

#' @title Funkcija za određivanje druge kontrolne veličine na osnovu jedne 
#' putanje geometrijskog Braunovog kretanja. 
#'
#' @param putanja numeric. Vektor koji predstavlja jednu putanju geometrijskog
#' Braunovog kretanja. Cijena evropske kol opcije zavisi samo od cijene akcije
#' u trenutku dospijeća, a ne od prosjeka svih cijena do tog trenutka.
#'
#' @return cijena evropske kol opcije dobijena simulacijom trajektorije procesa
#' kretanja cijena akcija
KontrolnaTip2 <- function(putanja) {
  exp(-r * T) * max(putanja[length(putanja)] - K, 0)
}

#' @title Blek-Šolsova formula
#' 
#' @description Funkcija koja implementira Blek-Šolsovu formulu za cijenu 
#' evropske kol opcije u trenutku t=0 sa ugovorenom cijenom K i vremenom do 
#' isteka opcije T. Njena povratna vrijednost je matematičko očekivanje dobiti 
#' evropske kol opcije pri valuaciji neutralnoj od rizika, odnosno matematičko 
#' očekivanje druge kontrolne promjenljive.
#' 
#' @param S_0 numeric.  Početna cijena akcije
#' @param K numeric. Ugovorena cijena opcije
#' @param T numeric. Ugovoreno vrijeme
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#'
#' @return B-Š cijena evropske opcije
#'
#' @examples
#' #cijena_bs <- AV_GBM_v(100, 0.05, 0.3, 1, 100, 1)
EvropskaOpcija <- function(S_0, K, T, r, sigma) {
  dl <- (log(S_0 / K) + (r + sigma ^ 2 / 2) * T) / (sqrt(T) * sigma)
  d2 <- dl - sqrt(T) * sigma
  cijenaEV <- S_0 * pnorm(dl) - K * exp(-r * T) * pnorm(d2) 
  return(cijenaEV)
}

E_KontrolnaTip2 <- EvropskaOpcija(S_0, K, T, r, sigma)

#' @title Funkcija za određivanje treće kontrolne slučajne veličine na osnovu 
#' putanje geometrijskog Braunovog kretanja. 
#' 
#' @description Primjenjuje se martingalski pristup po kom
#' je cijena opcije diskontovana očekivana vrijednost dobiti opcije pri
#' bezrizičnoj kamatnoj stopi.
#'
#' @param putanja Vektor koji predstavlja jednu putanju geometrijskog
#' Braunovog kretanja. Cijena geometrijske azijske opcije zavisi od geometrijske
#' sredine cijena akcija u intervalu [0,T]
#'
#' @return cijena geometrijske azijske opcije dobijena simulacijom trajektorije
#' procesa kretanja cijena akcija
KontrolnaTip3 <- function(putanja) {
  exp(-r * T) * max((prod(putanja[-1])) ^ (1 / n) - K, 0)
}

#' @title Cijena geometrijke azijske opcije pri mjeri neutralnoj od rizika
#' 
#' @description Funkcija koja implementira formulu za računanje cijene 
#' geometrijske azijske opcije čija cijena je data eksplicitnom formulom.
#' Posmatra se opcija u trenutku t=0. 
#' 
#' @param S_0 numeric. Početna cijena akcije
#' @param K numeric. Ugovorena cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka

#' @return cijena geometrijske azijske opcije
#'
#' @examples cijena_geom <- AzijskaOpcija_geom(100, 70, 0.05, 0.3, 100)
AzijskaOpcija_geom <- function(S_0, K, r, sigma, n) {
  # predefinisane konstante
  delta_t <- T / n
  c_1 <- r - sigma ^ 2 / 2
  a <- log(S_0) + c_1 * delta_t + (c_1 * (T - delta_t)) / 2
  b <- sigma ^ 2 * delta_t + sigma ^ 2 * (T - delta_t) * (2 * n - 1) / (6 * n)
  x <- (a - log(K) + b) / sqrt(b)
  cijena_geom <- exp(-r * T) * (exp(a + b / 2) * pnorm(x) - K * pnorm(x - sqrt(b)))
  return (cijena_geom)
}

# Matematičko očekivanje treće kontrolne slučajne veličine. 
E_KontrolnaTip3 <- AzijskaOpcija_geom(S_0, K, r, sigma, n)

# Prethodno definisana funkcija za generisanje geometrijskog Braunovog kretanja.
GB_v <- function(S_0, r, sigma, T, n, N) {
  delta_t <- T / n
  c_1 <- (r - (sigma ^ 2 / 2)) * delta_t
  c_2 <- sqrt(delta_t) * sigma
  norm <- matrix(rnorm(n * N, c_1, c_2), nrow = N, ncol = n)
  log_prirastaj <- cbind(log(S_0), norm)
  log_putanja <- t(apply(log_prirastaj, 1, cumsum))
  return(exp(log_putanja))
}

#funckija AzijskaOpcijaCV koristi CV metod pri racunanju cene azijske opcije
#' @title Metoda kontrolnih promjenljivih za određivanje cijene azijske opcije
#' 
#' @param S_0 numeric. Početna cijena akcije
#' @param K numeric. Ugovorena cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka
#' @param N numeric. Broj trajektorija geometrijskog Braunovog kretanja 
#' @param KontrolnaTip closure. Funkcija kojom je definisana željena
#' kontrolna slučajna veličina
#' @param E_KontrolnaTip numeric. Matematičko očekivanje odabrane kontrolne
#' slučajne veličine
#' @param br_kontrolnih_simulacija numeric. Broj pomoćnih simulacija za 
#' određivanje parametra c
#' @param beta numeric. Nivo povjerenja za konstrukciju intervala povjerenja
#'
#' @return cijena aritmetičke azijske opcije koristeći metodu kontrolnih
#' promjenljivih za redukciju disperzije
#'
#' @examples
#' cijena_azijska_aritmeticka <- AzijskaOpcijaCV(S_0, K, T, r, sigma, n, N, 
#' KontrolnaTip3, E_KontrolnaTip3, br_kontrolnih_simulacija, beta)
AzijskaOpcijaCV<-function(S_0, K, T, r, sigma, n, N, KontrolnaTip,
                          E_KontrolnaTip, br_kontrolnih_simulacija, beta){
  # dodatne ("pilot") simulacije za određivanje kontrolnog parametra c
  # vektor kontrolnih slučajnih veličina (Y)
  kontrolne_c <- rep(0, br_kontrolnih_simulacija)
  
  # vektor cijena aritmetičke azijske opcije dobijen kontrolnim simulacijama (X)
  cene_opcije_aritm <- rep(0, br_kontrolnih_simulacija)
  for (i in 1:br_kontrolnih_simulacija) {
    # jedna putanja geometrijskog Braunovog kretanja
    putanja <- GB_v(S_0, r, sigma, T , n, 1)
    # vrijednost kontrolne slučajne veličine za datu putanju
    kontrolne_c[i] <- KontrolnaTip(putanja)
    # cijena za datu putanju koristeći izraz za dobit opcije
    cene_opcije_aritm[i] <- exp(-r * T) * max(mean(putanja[-1]) - K, 0)
  }
  # primjenjujemo robusnu regresiju za određivanje kontrolnog parametra c
  # model: X = beta0 + beta1*Y + epsilon
  # napomena: zamijenjene su uloge X i Y u odnosu na standardnu notaciju 
  # linearnog modela
  # optimalno c: c_zvijezda = -beta1, gde je beta1 jedan od koeficijenata modela 
  # robusne regresije
  c_zvijezda<- -(rlm(cene_opcije_aritm~kontrolne_c)$coefficients[2])
  
  # Monte Karlo simulacije
  # vektor kontrolnih slučajnih veličina za MK metod (Y)
  kontrolne <- rep(0, N)
  # vektor cijena aritmetičke azijske opcije (X)
  cijene_opcije_CV <- rep(0, N)
  
  for (i in 1:N) {
    # jedna putanja geometrijskog Braunovog kretanja
    putanja <- GB_v(S_0, r, sigma, T , n, 1)
    # primjenjujemo odabranu funkciju KontrolnaTip za određivanje kontrolne
    # slučajne veličine za datu putanju
    kontrolne[i] <- KontrolnaTip(putanja)
    
    # cijena azijske opcije zavisi od prosjeka cijena u intervalu [0,T]
    cijena_az_opcije <- exp(-r * T) * max(0, mean(putanja[-1]) - K)
    
    # ocjena za cijenu aritmetičke azijke opcije koristeći metodu kontrolnih 
    # promjenljivih
    cijene_opcije_CV[i] <-
      cijena_az_opcije + c_zvijezda * (kontrolne[i] - E_KontrolnaTip)
  }
  # standardna devijacija niza ocjena
  std <- sd(cijene_opcije_CV)
  
  # konačna MK ocjena je prosjek dobijenih cijena
  CijenaKonacno_MK <- mean(cijene_opcije_CV) 
  StudentovTest <- t.test(cijene_opcije_CV, conf.level = beta)
  IP <- StudentovTest$conf.int
  greska <-
    100 * (1 / (as.numeric(StudentovTest$estimate))) *
    ((StudentovTest$conf.int[2] - StudentovTest$conf.int[1]) / 2)
  cat("Tackasta ocjena za cijenu aritmeticke azijske opcije koriscenjem metode kontrolnih promjenljivih je: ",
      CijenaKonacno_MK, "\n")
  cat("Standardna devijacija je: ", std, "\n")
  cat(beta * 100, "% interval poverenja je: ", IP, "\n")
  cat("Greska ocjene u procentima je: ", greska, "\n")
}
AzijskaOpcijaCV(S_0, K, T, r, sigma, n, N, KontrolnaTip3, E_KontrolnaTip3,
                br_kontrolnih_simulacija, beta)
