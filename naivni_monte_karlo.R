# funkcija kojom se generiše geometrijsko Braunovo kretanje
# preuzeta iz fajla "GBM.R"
GBM_v <- function(S_0, r, sigma, T, n, N) {
  # podjela vremenskog intervala
  delta_t <- T / n
  
  # izračunavamo unaprijed poznate konstante
  c_1 <- (r - (sigma ^ 2 / 2)) * delta_t
  c_2 <- sqrt(delta_t) * sigma
  
  # matrica u kojoj smještamo slučajne veličine iz normalne raspodjele
  norm <- matrix(rnorm(n * N, c_1, c_2), nrow = N, ncol = n)
  
  # kreiramo novu matricu u koju je dodata kolona početnih vrijednosti cijene
  log_prirastaj <- cbind(log(S_0), norm)
  
  # računamo kumulativnu sumu vrijednosti matrice log_prirastaj za svaku vrstu
  # kako funkcija apply za argument MARGIN=1 računa kumulativnu sumu svake vrste
  # po kolonama i smješta vrijednosti u matricu log_putanja po vrstama obrnu se
  # dimenzije početne matrice
  # potrebno je transponovati dobijenu matricu kako bismo dobili matricu u kojoj 
  # svaka vrsta predstavlja jednu trajektoriju geometrijskog Braunovog kretanja
  log_putanja <- t(apply(log_prirastaj, 1, cumsum))
  return (exp(log_putanja))
}

#' @title Funkcija koja koristi obični (naivni) Monte-Karlo metod za 
#' računanje cijene aritmetičke azijske opcije.
#' 
#' @description Za vrednovanje opcije koristi se princip valuacije neutralan od 
#' rizika. Po ovom principu cijena opcije je diskontovana očekivana vrijednost 
#' dobiti opcije pri bezrizičnoj kamatnoj stopi.
#' Monte-Karlo ocjena za cijenu je prosjek cijena opcije dobijenih 
#' u simulacijama.
#' 
#' @param K numeric. Ugovorena cijena akcije na koju se opcija odnosi
#' @param S_0 numeric. Početna cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka
#' @param N numeric. Broj Monte-Karlo simulacija (ujedno i broj generisanih 
#' trajektorija)
#' @param beta numeric. Nivo povjerenja za konstrukciju intervala povjerenja
#'
#' @examples
#' AzijskaOpcija(70, 100, 0.05, 0.3, 1, 100, 10000, 0.95)

AzijskaOpcija <- function(K, S_0, r, sigma, T, n, N, beta) {
  cijena_opcije <- 0
  
  # vektor u koji smještamo cijene opcije dobijene martingalskim pristupom
  # računa se cijena opcije za svaku simulaciju, odnosno trajektoriju
  cijene <- rep(0, N)
  
  # faktor diskontovanja je konstantna vrijednost
  faktor_diskontovanja <- exp(-r * T)
  
  for (i in 1:N) {
    # jedna trajektorija geometrijskog Braunovog kretanja
    putanja <- GBM_v(S_0, r, sigma, T, n, 1)
    # cijena akcije se računa po sledećoj formuli
    cijene[i] <- faktor_diskontovanja * max(0, mean(putanja[-1]) - K)
  }
  # standardna devijacija niza dobijenih cijena opcije
  std <- sd(cijene)
  # Monte-Karlo ocjena cijene opcije
  cijena_opcije <- mean(cijene)
  # primjenjujemo studentov test na niz cijena da bismo dobili 
  # interval povjerenja
  studentov_test <- t.test(cijene, conf.level = beta)
  interval_povjerenja <- studentov_test$conf.int
  # računamo grešku Monte-Karlo ocjene
  greska <-
    100 * (1 / (as.numeric(studentov_test$estimate))) *
    ((studentov_test$conf.int[2] - studentov_test$conf.int[1]) / 2)
  # štampamo vrijednost ocjene, standardnu devijaciju, interval povjerenja i 
  # samu grešku
  cat("Tačkasta ocena za cenu azijske opcije je: ", cijena_opcije, "\n")
  cat("Standardna devijacija je: ", std, "\n")
  cat(beta * 100, "% interval povjerenja je: ", interval_povjerenja, "\n")
  cat("Greška ocjene izražena u procentima je: ", greska, "\n")
}
# primjer
set.seed(77)
AzijskaOpcija(100, 100, 0.05, 0.4, 1, 12, 1000, 0.95)
# ako nas zanima vrijeme potrebno da bi se kod izvrsio u sekundama)
# system.time(AzijskaOpcija(30, 50, 0.05, 0.3, 1, 100, 10000, 0.95))

