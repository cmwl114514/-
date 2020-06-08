result<-main(step = 90000,burnin = 10000,
             N=3*10^3,
             ProbA=c(0.7,0.7),
             ProbB=c(0.6,0.6),
             ProbC=c(0.3,0.3),
             ProbD=c(0.5,0.5),
             ProbX=0.6)

var_list<-lapply(result, var)
E_list<-lapply(result, mean)
M_list<-lapply(result, median)