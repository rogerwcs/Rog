bsolver1=function(N,lambda,b0=1,bN=1)
{
  det=(1-exp(-lambda))^0.5
  u1=(1+det)
  u2=(1-det)
  sol1=u1^c(0:N)
  sol2=u2^c(0:N)
  if (lambda==0) sol2=sol2*c(0:N)
  # c1 sol1[1]+c2 sol2[1]=b0
  #c1 sol1[N+1]+c2 sol2[N+1]=bN
  A=matrix(c(sol1[1],sol2[1],sol1[N+1],sol2[N+1]),nrow=2,ncol=2, byrow=TRUE)
  cc=solve(A)%*%c(b0,bN)
  b=cc[1]*sol1+cc[2]*sol2
  return(b[2:N])
  %teste
}

#This calculates moments up to order two of \tau(1)

#When b0=bN=1: counting number steps by "1"until absobed or absorbs

#When b0=1,bN=0: on the event that "1" was not absorbed
#When b0=0,bN=1: on the event that "1" did not absorb
# In either case the zeroth moment is the probability of the event
# and the other moments are the expected number of steps on this event

Moments=function(N,b0=1,bN=1)
{
  f=function(x) return(bsolver(N,x,b0,bN)[1])
  e0=f(0)
  e1=grad(f,0,side=c(1))
  g=function(x) {if (x==0) return(0)
    else return((f(x)-e0-e1*x)*2/x)}
  
  e2=grad(g,0,side=c(1))
  return(c(e0,-e1,e2))	
}
Ns=c(3:15)
#The second row in each of E11 and E10 is the expectation
#The third row is the second moment

E11=sapply(Ns,function(x) Moments(x))
E10=sapply(Ns,function(x) Moments(x,1,0))
