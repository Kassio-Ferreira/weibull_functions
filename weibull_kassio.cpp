#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec powvec(vec vetor, float expoente){
  vec saida(vetor.n_elem);
  for(int i=0; i<vetor.n_elem; i++){
    saida(i) = pow(vetor(i),expoente);
  }
  return saida;
}

// [[Rcpp::export]]
vec dweibull_standard(vec t, float alpha, float beta) {
  return beta*powvec(t,beta-1)/pow(alpha,beta)%exp(-powvec(t/alpha,beta));
}

// [[Rcpp::export]]
vec weibull_survivor(vec t, float alpha, float beta){
  return exp(-powvec(t/alpha,beta));
}

// [[Rcpp::export]]
vec pweibullC(vec t, float alpha, float beta){
  return 1 - weibull_survivor(t, alpha, beta);
}

// [[Rcpp::export]]
vec hazard_weibull(vec t, float alpha, float beta){
  return -log(1-pweibullC(t,alpha,beta));
}

// [[Rcpp::export]]
float fgm_weibull(int s, float alpha, float beta){
  return pow(alpha,s)*exp(lgamma(1 + s/beta));
}

// gerador da distribuição weibull de dois parâmetros pelo método da inversa
// [[Rcpp::export]]
vec rweibullC(float n, float beta, float alpha){
  return pow(log(1/randu(n)),(1/beta)) * alpha;
}

// gerador da distribuição weibull de dois parâmetros pelo método da aceitação-rejeicao
// [[Rcpp::export]]
vec rweibullAR(int n, float beta, float alpha, float a, float b, float k){
  vec amostraFinal(n);
  int cont = 0;
  float X, Y, fx;
  
  while(cont < n){
    X = a + (b-a)*conv_to<float>::from(randu(1));
    Y = k*conv_to<float>::from(randu(1));
    fx = beta*pow(X,beta-1)/pow(alpha,beta)*exp(-pow(X/alpha,beta));
      
    if(Y < fx){
      amostraFinal(cont) = X;
      cont++;
    }
    
  }
  return amostraFinal; 
}


/*** R
# densidade weibull
x = seq(0,3,0.05)
plot(x,dweibull_standard(x,2,2), t="l", lwd=6, col="red")
lines(x,dweibull(x,2,2), col="blue")

# Acumulada weibull
plot(x,pweibullC(x,2,2), lwd=6, col="red", t="l")
lines(x,pweibull(x,2,2), lwd=1, col="blue")

# hazard weibull
plot(x,hazard_weibull(x,2,2), lwd=6, col="red", t="l")
lines(x, -log(1-pweibull(x,2,2)), col="blue", t="l")

# função geradora de momentos
fgm_weibull(1,2,2)

# Gerador aleatório da Weibull padrão (método da inversa):
rweibullR = function(n,beta,alpha) log(1/runif(n,0,1))^(1/beta)*alpha #implementado em R

par(mfrow=c(1,3))
hist(rweibullR(100,6,2))
hist(rweibull(100,6,2))
hist(rweibullC(100,6,2))

# comparação de amostras geradas pela função do C e R-base
N = 100
alpha=2; beta=6
cont = 0
for(i in 1:N){
  if(ks.test(rweibull(100,beta,alpha), rweibullC(100,beta,alpha))$p.value < 0.05)
    cont = cont + 1
}

cont/N

# tempo de processamento
d = microbenchmark::microbenchmark(C=rweibullC(1000,4,5), R=rweibullR(1000,4,5), R_base = rweibull(1000,4,5))
print(d)

# Gerador aleatório da Weibull padrão (método da aceitação-rejeição):
geradorRejeicao <- function(N,beta,alpha,a,b,k){
  amostraFinal = NULL;

  while(length(amostraFinal) < N){

    X = runif(1,a,b)
    Y = runif(1,0,k)
    fx = dweibull_standard(X, alpha, beta)
  
    if(Y<fx) amostraFinal = c(amostraFinal,X)
  }
  return(amostraFinal);
}


par(mfrow=c(1,3))
alpha = 4; beta = 3
hist(rweibull(1000,beta,alpha), main="R-base", breaks=30)
hist(geradorRejeicao(1000,beta,alpha,0,8,3), main="R", breaks=30)
hist(rweibullAR(1000,beta,alpha,0,8,3), main="C", breaks=30)

# comparação de amostras geradas pela função do R e R-base
cont = 0
for(i in 1:N){
  if(ks.test(rweibull(100,beta,alpha), geradorRejeicao(100,beta,alpha,0,8,k=3))$p.value < 0.05)
    cont = cont + 1
}

cont/N

# comparação de amostras geradas pela função do C e R-base
cont = 0
for(i in 1:N){
  if(ks.test(rweibull(1000,alpha,beta), rweibullAR(100,alpha,beta,0,8,3))$p.value < 0.05)
    cont = cont + 1
}

cont/N

*/
