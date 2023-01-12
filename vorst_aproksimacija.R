# Vorstova aproksimacija za cijenu aritmetičke azijske opcije

# C_G - cijena geometrijske azijske opcije 
# C_A - cijena aritmetičke azijske opcije
# C_V - Vorstova aproksimacija za cijenu aritmetičke azijske opcije
# donja_granica - donja granica za cijenu aritmetičke azijske opcije
# gornja_granica - gornja granica za cijenu aritmetičke azijske opcije

#' @title Očekivanje prosjeka cijena akcija pri valuaciji neutralnoj od rizika
#'
#' @param S_0 numeric. Početna cijena akcije
#' @param K numeric. Ugovorena cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param T numeric. Vrijeme do isteka opcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo slučajni
#' proces, odnosno broj koraka
#'
#' @return traženo matematiko očekivanje
E_A_T <- function(S_0, K, r, sigma, T, n) {
  delta_t <- T / n
  return((S_0 / n) * exp(r * delta_t) * (1 - exp(r * delta_t * n)) / (1 - exp(r * delta_t)))
}

E_G_T <- function(S_0, K, r, sigma, T, n) {
  # predefinisne konstante
  c_1 <- r - sigma ^ 2 / 2
  delta_t <- T / n
  
  # a = E(ln(G(T)))
  a <- log(S_0) + c_1 * delta_t + (c_1 * (T - delta_t)) / 2
  
  # b = D(ln(G(T)))
  b <- sigma ^ 2 * delta_t + sigma ^ 2 * (T - delta_t) * (2 * n - 1) / (6 * n)
  
  # E(G(T)) u odnosu na mjeru neutralnu od rizika
  return(exp(a + b / 2))
}

#' @title Funkcija za određivanje cijene geometrijske azijske opcije po 
#' poznatoj formuli. 
#' 
#' @description Formula je Blek-Šolsovog tipa zasnovana na činjenici da je
#' proizvod slučajnih veličina sa log-normalnom raspodjelom slučajna veličina
#' koja takođe ima log-normalnu raspodjelu. Povratna vrijednost ove funkcije
#' je ujedno i donja granica za cijenu odgovarajuće aritemtičke azijske opcije.
#' 
#' @param S_0 numeric. Početna cijena akcije
#' @param K numeric. Ugovorena cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param T numeric. Vrijeme do isteka opcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo slučajni
#' proces, odnosno broj koraka
#'
#' @return cijena geometrijske azijske opcije
AzijskaOpcija_geom <- function(S_0, K, r, sigma, T, n) {
  # predefinisane konstante
  delta_t <- T / n
  c_1 <- r - sigma ^ 2 / 2
  a <- log(S_0) + c_1 * delta_t + (c_1 * (T - delta_t)) / 2
  b <- sigma ^ 2 * delta_t + sigma ^ 2 * (T-delta_t) * (2 * n - 1) / 6 / n
  x <- (a - log(K) + b) / sqrt(b)
  
  # formula za cijenu geometrijske azijske opcije u trenutku t = 0
  cijena_geom <-
    exp(-r * T) * (exp(a + b / 2) * pnorm(x) - K * pnorm(x - sqrt(b)))
  return (cijena_geom)
}

# funkcija kojom se računa gornja granica za cijenu aritemtičke azijske opcije
# isti parametri kao u prethodnoj funkciji
gornja_granica <- function(S_0, K, r, sigma, T, n) {
  E_A <- E_A_T(S_0, K, r, sigma, T, n)
  E_G <- E_G_T(S_0, K, r, sigma, T, n)
  cijena_geom <- AzijskaOpcija_geom(S_0, K, r, sigma, T, n)
  # formula za gornju granicu za cijenu aritmetičke azijske opcije
  return(exp(-r * T) * (E_A - E_G) + cijena_geom)
}

# Vorstova aproksimacija se dobija primjenom formule za cijenu geometrijske 
# azijske opcije sa korigovanom ugovorenom cijenom K'=K-[E(A(T))-E(G(T))]
# očekivanje se računa u odnosu na mjeru neutralnu od rizika
vorst_aproksimacija <-
  function(S_0, K, r, sigma, T, n, AzijskaOpcija_geom, gornja_granica) {
    # korigovana ugovorena cijena
    K1 <- K - (E_A_T(S_0, K, r, sigma, T, n) - E_G_T(S_0, K, r, sigma, T, n))
    
    # vorstova aproksimacija za cijenu aritmetičke azijske opcije
    vorst_cijena <- AzijskaOpcija_geom(S_0, K = K1, r, sigma, T, n)
    
    # ispisujemo vrijednosti Vorstove aproksimacije, kao i odgovarajuću gornju
    # i donju granicu
    cat("Vorstova aproksimacija za cijenu aritmetičke azijske opcije je: ",
        vorst_cijena, "\n")
    cat("Gornja granica za cijenu aritmetičke azijske opcije je: ",
        gornja_granica(S_0, K, r, sigma, T, n), "\n")
    cat("Donja granica za cijenu aritmetičke azijske opcije je: ", 
        AzijskaOpcija_geom(S_0, K = K, r, sigma, T, n), "\n")
  }

vorst_aproksimacija(100, 100, 0.1, 0.5, 1, 100, 
                    AzijskaOpcija_geom, gornja_granica)
