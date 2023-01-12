# Linearna kombinacija triju kontrolnih slučajnih veličina Y_1, Y_2, Y_2.
# Ocjenjujemo nepoznati parametar teta = E(X), X korelisana sa Y_1, Y_2, Y_3.
# Posmatramo linearni model X = beta_0 + beta_1*Y_1 + beta_2*Y2 + beta_3*Y_3

# Slučajne veličine Y_i, i=1,2,3 su redom:
# 1) Aritmetičkа sredina cijena akcija.
# 2) Evropska kol opcija čija je cijena eksplicitno data B-Š formulom.
# 3) Geometrijska azijska opcija za čiju cijenu postoji eksplicitna formula.
# Potrebno je navesti vrijednosti parametara na početku, jer je prije 
# pokretanja glavne funkcije AzijskaOpcija_CV_3 potrebno izračunati matematička
# očekivanja kontrolnih slučajnih veličina za zadate parametre.

library(MASS) 
S_0 <- 50
K <- 50
r <- 0.05
sigma <- 0.6
T <- 1
n <- 12
N <- 5000
delta_t <- T / n
beta <- 0.95
br_kontrolnih_simulacija <- 100

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
# E_KontrolnaTip2 - matematičko očekivanje druge kontrolne slučajne veličine
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
#' Posmatra se opcija u trenutku t = 0. 
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
  # predefinisane konstante delta_t, c_1, a, b, x
  delta_t <- T / n
  c_1 <- r - sigma ^ 2 / 2
  a <- log(S_0) + c_1 * delta_t + (c_1 * (T - delta_t)) / 2
  b <- sigma ^ 2 * delta_t + sigma ^ 2 * (T - delta_t) * (2 * n - 1) / (6 * n)
  x <- (a - log(K) + b) / sqrt(b)
  
  # formula za cijenu geometrijske azijske opcije
  cijena_geom <-
    exp(-r * T) * (exp(a + b / 2) * pnorm(x) - K * pnorm(x - sqrt(b)))
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

#' @title Uopštena metoda kontrolnih promjenljivih za određivanje cijene 
#' aritmetičke azijske opcije
#' 
#' @description Kombinuju se tri kontrolne slučajne veličine: 
#' aritmetička sredina cijena akcija, evropska kol opcija i geometrijska 
#' azijska opcija. 
#'
#' @param S_0 numeric. Početna cijena akcije na koju se opcija odnosi
#' @param K numeric. Ugovorena cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka
#' @param N numeric. Broj trajektorija geometrijskog Braunovog kretanja 
#' @param E_KontrolnaTip1 numeric. Matematičko očekivanje prve kontrolne
#' slučajne veličine
#' @param E_KontrolnaTip2 numeric. Matematičko očekivanje druge kontrolne
#' slučajne veličine
#' @param E_KontrolnaTip3 numeric. Matematičko očekivanje treće kontolne 
#' slučajne veličine. 
#' @param br_kontrolnih_simulacija numeric. Broj pomoćnih simulacija za 
#' određivanje parametra c
#' @param beta numeric. Nivo povjerenja za konstrukciju intervala povjerenja
#'
#' @return cijena aritmetičke azijske opcije koristeći uopštenu metodu 
#' kontrolnih promjenljivih za redukciju disperzije
#'
#' @examples
#' cijena_az_cv3 <- AzijskaOpcijaCV(S_0, K, T, r, sigma, n, N, E_KontrolnaTip1, 
#'                                  E_KontrolnaTip2, E_KontrolnaTip3, 
#'                                  br_kontrolnih_simulacija, beta)

AzijskaOpcijaCV_3 <- function(S_0, K, T, r, sigma, n, N,
                              E_KontrolnaTip1, E_KontrolnaTip2, E_KontrolnaTip3,
                              br_kontrolnih_simulacija, beta){
  # dodatne ("pilot") simulacije za određivanje kontrolnog parametra c
  # u opštem slučaju imamo vektor sa tri kontrolna parametra
  
  # vektori kontrolnih slučajnih veličina (Y)
  kontrolne_c_1 <- rep(0, br_kontrolnih_simulacija)
  kontrolne_c_2 <- rep(0, br_kontrolnih_simulacija)
  kontrolne_c_3 <- rep(0, br_kontrolnih_simulacija)
  
  # vektor cijena aritmetičke azijske opcije dobijen kontrolnim simulacijama (X)
  cijene_opcije_aritm <- rep(0, br_kontrolnih_simulacija)
  for (i in 1:br_kontrolnih_simulacija) {
    # jedna putanja geometrijskog Braunovog kretanja
    putanja <- GB_v(S_0, r, sigma, T , n, 1)
    # vrijednosti kontrolnih slučajnih veličina za datu putanju
    kontrolne_c_1[i] <- KontrolnaTip1(putanja)
    kontrolne_c_2[i] <- KontrolnaTip2(putanja)
    kontrolne_c_3[i] <- KontrolnaTip3(putanja)
    # cijena za datu putanju koristeći izraz za dobit opcije
    cijene_opcije_aritm[i] <- exp(-r * T) * max(mean(putanja[-1]) - K, 0)
  }
  # primjenjujemo robusnu regresiju za određivanje optimalnog parametra
  # c_zvijezda = (c_1, c_2, c_3)
  # model: X = beta0 + beta1*Y_1 + beta2*Y_2 + beta3*Y_3 + epsilon
  # napomena: zamijenjene su uloge X i Y u odnosu na standardnu notaciju 
  # linearnog modela
  # optimalno c: c_i_zvijezda = -beta_i, i = 1, 2, 3, gde je beta_i 
  # odgovarajući koeficijent modela robusne regresije
  model_rlm <- rlm(cijene_opcije_aritm~kontrolne_c_1+kontrolne_c_2+kontrolne_c_3,
                   maxit = 50)
  c_zvijezda_1<- -model_rlm$coefficients[2]
  c_zvijezda_2<- -model_rlm$coefficients[3]
  c_zvijezda_3<- -model_rlm$coefficients[4]

  # Monte Karlo simulacije
  # vektori kontrolnih slučajnih veličina za MK metod (Y)
  kontrolne1 <- rep(0, N)
  kontrolne2 <- rep(0, N)
  kontrolne3 <- rep(0, N)
  # vektor cijena aritmetičke azijske opcije (X)
  cijene_opcije_CV <- rep(0, N)
  
  for (i in 1:N) {
    # jedna putanja geometrijskog Braunovog kretanja
    putanja <- GB_v(S_0, r, sigma, T , n, 1)
    # realizacije slučajnih veličina za datu putanju
    kontrolne1[i] <- KontrolnaTip1(putanja)
    kontrolne2[i] <- KontrolnaTip2(putanja)
    kontrolne3[i] <- KontrolnaTip3(putanja)
    
    # cijena azijske opcije zavisi od prosjeka cijena u intervalu [0,T]
    cijena_az_opcije <- exp(-r * T) * max(0, mean(putanja[-1]) - K)
    
    # ocjena za cijenu aritmetičke azijke opcije koristeći uopštenu metodu 
    # kontrolnih promjenljivih
    cijene_opcije_CV[i] <-
      cijena_az_opcije + c_zvijezda_1 * (kontrolne1[i] - E_KontrolnaTip1) +
      c_zvijezda_2 * (kontrolne2[i] - E_KontrolnaTip2) + c_zvijezda_3 * 
      (kontrolne3[i] - E_KontrolnaTip3)
  }
  # konačna MK ocjena
  CijenaKonacno_MK <- mean(cijene_opcije_CV) 
  # standardna devijacija niza ocjena
  std <- sd(cijene_opcije_CV)
  StudentovTest <- t.test(cijene_opcije_CV, conf.level = beta)
  IP <- StudentovTest$conf.int
  greska <-
    100 * (1 / (as.numeric(StudentovTest$estimate))) *
    ((StudentovTest$conf.int[2] - StudentovTest$conf.int[1]) / 2)
  cat("Uopstena metoda kontrolnih promjenljivih za redukciju disperzije \n")
  cat("Tackasta ocjena za cijenu aritmeticke azijske opcije je: ",
      CijenaKonacno_MK, "\n")
  cat("Standardna devijacija je: ", std, "\n")
  cat(beta * 100, "% interval poverenja je: ", IP, "\n")
  cat("Greska ocjene u procentima je: ", greska, "\n")
}
set.seed(77)
AzijskaOpcijaCV_3(S_0, K, T, r, sigma, n, N, E_KontrolnaTip1, E_KontrolnaTip2, 
                  E_KontrolnaTip3, br_kontrolnih_simulacija, beta)
