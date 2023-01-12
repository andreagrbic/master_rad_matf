#' @title Funkcija za generisanje geometrijskog Braunovog kretanja koja je 
#' prilagođena metodi antitetičkog uzorkovanja za redukciju disperzije 
#' (Antithetic Variates). Generišu se putanje na odgovarajući način.
#' 
#' @description Izbjegava se korišćenje ugnježđenih for petlji,
#' funkcija je vektorizovana. 
#' 
#' Generisanje uzorka iz normalne raspodjele se dijeli na dva dijela.
#' Ideja je da se generiše uzorak tako da je prva polovina uzorka negativno 
#' korelisana sa drugom polovinom i da su obje polovine uzorka iz normalne 
#' N(c_1, c_2) raspodjele. Generisani parovi trajektorija izgledaju kao slike u
#' ogledalu.
#' 
#' @param S_0 numeric. Početna cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka
#' @param N numeric. Broj trajektorija geometrijskog Braunovog kretanja 
#'
#' @return Matrica trajektorija
#'
#' @examples
#' AV_GBM_v(100, 0.05, 0.3, 1, 100, 1)
AV_GBM_v <- function(S_0, r, sigma, T, n, N) {
  # podjela vremenskog intervala
  delta_t <- T / n
  
  # izračunavamo unaprijed poznate konstante
  c_1 <- (r - (sigma ^ 2 / 2)) * delta_t
  c_2 <- sqrt(delta_t) * sigma
  
  # u matricu Z sa N/2 vrsta smještamo slučajne veličine iz normalne raspodjele
  Z <- matrix(rnorm(n * ceiling(N / 2), c_1, c_2),
              nrow = ceiling(N / 2), ncol = n)
  
  # u matrici av_norm na matricu Z nadovezujemo odgovarajući antitetički uzorak
  # matrica av_norm sadrži slučajne veličine Z i 2*c_1-Z što su traženi 
  # antitetički parovi replikacija slučajnih veličina iz N(c_1, c_2) raspodjele
  # slučajne veličine unutar svakog para su negativno korelisane
  av_norm <- rbind(Z, 2 * c_1 - Z)
  
  # u matricu av_log_prirastaj je dodata kolona početnih vrijednosti 
  # cijene akcije
  av_log_prirastaj <- cbind(log(S_0), av_norm)
  
  # računamo kumulativnu sumu vrijednosti matrice log_prirastaj za svaku vrstu
  # kako funkcija apply za argument MARGIN=1 računa kumulativnu sumu svake vrste
  # po kolonama i smješta vrijednosti u matricu log_putanja po vrstama obrnu se
  # dimenzije početne matrice
  # potrebno je transponovati dobijenu matricu kako bismo dobili matricu u kojoj 
  # svaka vrsta predstavlja jednu trajektoriju geometrijskog Braunovog kretanja
  log_putanja <- t(apply(av_log_prirastaj, 1, cumsum))
  return(exp(log_putanja))
}

# izgled trajektorija sa funkcijom AV_GBM_v
# ovim kodom ilustrujemo kako izgledaju trajektorije geometrijskog Braunovog 
# kretanja kada koristimo replikacije slučajnih veličina iz normalne raspodjele
# koje su negativno korelisane


set.seed(777)
putanje <- AV_GBM_v(50, 0.1, 0.2, 1, 300, 15)
boje <- rainbow(14)
x <- c(1:301)
y_granice <- c(min(putanje), max(putanje))
plot(x, putanje[1, ], ylim = y_granice, type = 'l', 
     xlab = expression(paste(sigma," = 0.2")), ylab = "", lty = 1)
for (i in 2:15)
  lines(x, putanje[i, ], col = boje[i])

