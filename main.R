main<-function(loop=100,step=10^4,
               burnin=10^3,
               N=3*10^3,
               ProbA=c(0.3,0.7),
               ProbB=c(0.7,0.4),
               ProbC=c(0.5,0.3),
               ProbD=c(0.3,0.9),
               ProbX=0.6){
  s=Sys.time()
  N_list_EM<-simulation_EM(loop = loop,step = 60,N=N,ProbA=ProbA,ProbB=ProbB,ProbC=ProbC,ProbD=ProbD,ProbX=ProbX)
  N_list_mis<-simulation_MCMC_mis(loop = loop,step = step,burnin = burnin,N=N,ProbA=ProbA,ProbB=ProbB,ProbC=ProbC,ProbD=ProbD,ProbX=ProbX)
  N_list_com<-simulation_MCMC_complete(loop = loop,step = step,burnin = burnin,N=N,ProbA=ProbA,ProbB=ProbB,ProbC=ProbC,ProbD=ProbD,ProbX=ProbX)
  e=Sys.time()
  print(e-s)
  boxplot(N_list_EM,N_list_mis,N_list_com,names = c("EM","MCMC for missing data","MCMC for complete data"))
  abline(h=N,lwd=1,col="red",lty=2)
  return(list(N_list_EM,N_list_mis,N_list_com))
}

result<-main(step = 90000,burnin = 10000)
var_list<-lapply(result, var)
E_list<-lapply(result, mean)
M_list<-lapply(result, median)