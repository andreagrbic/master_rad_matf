#' Funkcija kojom se generiše geometrijsko Braunovo kretanje.
#' @description Izbjegava se korišćenje ugnježđenih for petlji,
#' funkcija je vektorizovana. 
#' 
#' Generišu se logaritmi priraštaja, odnosno razlike logaritama cijena akcija.
#' 
#' Povratna vrijednost funkcije je matrica čija svaka vrsta predstavlja
#' jednu trajektoriju geometrijskog Braunovog kretanja.
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
#' putanja <- GBM_v(100,0.05,0.3,1,100,1)
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

# primjer generisanja putanja i grafička reprezentacija
set.seed(7)
putanje1 <- GBM_v(50, 0.1, 0.1, 1, 300, 10)
putanje2 <- GBM_v(50, 0.1, 0.5, 1, 300, 10)
boje <- rainbow(9)
par(mfrow = c(1, 2))
x <- c(1:301)
y_granice <- c(min(min(putanje1), min(putanje2)),
               max(max(putanje1), max(putanje2)))
plot(x, putanje1[1, ], ylim = y_granice, type = 'l', 
     xlab = "sigma = 0.1", ylab = "", lty = 1)
for (i in 2:10)
  lines(x, putanje1[i, ], col = boje[i])
plot(x, putanje2[1, ], ylim = y_granice, type = 'l', 
     xlab = "sigma = 0.5", ylab = "", lty = 1)
for (i in 2:10)
  lines(x, putanje2[i, ], col = boje[i])
mtext("Trajektorije geometrijskog Braunovog kretanja sa volatilnošću 10% i 50%",
      side = 3, line = -2, outer = TRUE)
