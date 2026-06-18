import numpy as np
from models.campsb import CAMPSB, CAMPSB_var
from models.other import CTJAC, DCTJAC
import copy
from matplotlib import pyplot as plt

def rk45f(tk,xk,hk,hmin,hmax,tol,vfield,params): #vfield és DF? o què és en el nostre cas?
    '''
    Runge-Kutta Fehlberg orders 4-5.
    Input: tk (time), xk (state), hk (suggested time step forward (+) or
           backward (-), hmin (min allowed step), hmax (max allowed step),
           tol (truncation error), vfield (vectorfield, autonomous in our case)
    Output: tk1, xk1: time and state forward or backward in the propagation
            hk1: suggested time step for the next call
            err: truncation error found in the propagation from tk to tk1
    Note. It accepts a call with abs(hk)<hmin. In this case it proceeds to
          integrate with step h regardless of error.
    '''
    xk = np.array(xk)
    shp = xk.shape
    if len(shp)==1:
        pass
    else:
        nr,nc=xk.shape
        #if (nr<nc), xk=xk' end input xk can be either row or column vector
    error=10e10 #maxreal value, or a large value
    hs=hk
    if (abs(hk)>hmin):
        iz=0
    else:
        iz=1

    while (error > tol and iz<2):
        if (iz==1):#if iz=1 this will be last iteration (min step)
            iz=2
        ha=hs
        k1=ha*vfield(xk,params)
        k2=ha*vfield(xk+0.25*k1,params)
        k3=ha*vfield(xk+0.09375*k1+0.28125*k2,params)
        k4=ha*vfield(xk+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3,params)
        k5=ha*vfield(xk+(439/216)*k1-8*k2+(3680/513)*k3-(845/4104)*k4,params)
        k6=ha*vfield(xk-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-0.275*k5,params)
        error=np.linalg.norm(k1/360-(128/4275)*k3-(2197/75240)*k4+k5/50+(2/55)*k6)
        hp=abs(ha/hmax)
        if error<=tol*hp**5:
            hs=hmax*np.sign(ha)
        else:
            hs=0.9*ha*(tol/error)**0.2
            if (abs(hs)<=hmin):
                hs=hmin*np.sign(ha)
                iz=iz+1
    err=error
    hk1=hs
    tk1=tk+ha
    xk1=xk+(16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6
    #if (nr<nc), xk1=xk1' end
    return [tk1,xk1,hk1,err]

def cjdcj(xvec,params):
    '''
    Returns the Jacobi constant and its gradient (#TODO) given an initial point
    '''
    print("cjdcj")
    cjac = CJAC(xvec=xvec[:3],pvec=xvec[3:],params=params)
    dcjac = DCTJAC(xvec=xvec[:3],pvec=xvec[3:],params=params)
    return cjac, dcjac

def propTiTf(ti,xi,tf,vfield,h,hmin,hmax,tol,params):
    '''
    It propagates the i.c. (ti,xi) in the vectorfiled vfield till time tf
    and returns the final value xf (and the trajectory, if wanted).

    The vectorfield can include or not variational equations. If variational
    equations are included xi must also contain the initial condition for
    these equations (this is, xi is always an i.c. for vfield).

    h,hmin,hmax,tol are the usual parameters for the RKF integrator.
    Doesn't matter the sign of the initial step h. Integration always goes
    from ti to tf determining this way the initial (pos or neg) step.

    Returns:
        xfa: array of positions of the orbit
        tfa: array of times corresponding to each xfi from xfa
    '''
    TOLT = 1.e-12 #minimum time to integrate tolerance
    assert (abs(ti-tf) > TOLT), "propTiTf: initial time and final time are too close for propagation"
    h=abs(h)
    if (tf < ti):#enrere en el temps
        h = -h
    t=ti
    x=xi
    tfa=[]
    xfa=[]
    while ((t<tf and h>0) or (t>tf and h<0)):
        tfa.append(t)
        xfa.append(x)
        t,x,h,_ = rk45f(t,x,h,hmin,hmax,tol,vfield,params)
    while (abs(t-tf) > 1.e-12): #refinament de la solució per fer l'últim passet
        h=tf-t
        [t,x,_,_]=rk45f(t,x,h,hmin,hmax,tol,vfield,params)
    tfa.append(t)
    xfa.append(x)
    xfa = np.array(xfa) #cal agafar l'últim només no?
    tfa = np.array(tfa) #ídem
    #xfa = xfa.transpose()
    #plt.scatter(xfa[0],xfa[1])
    #plt.show()

    #print("xfa",xfa.shape,"\n",xfa)
    #print("tfa",tfa.shape,"\n",tfa)
    return [tfa, xfa]

def fdfArcFromTiToTf(ti,xi,tf,h,hmin,hmax,tol,vfield,hnd,params):
    '''
    It propagates the state (ti,xi) in the vectorfiled vfield till time tf
    and returns the final state xf and the state transition matrix stm
    computed either from variational equations or numerically.

    - If the variational equations are included in the vectorfield give hnd=0
      and the stm is obtained integrating them. Always, input xi should be
      the state, this is, without the i.c. of the variational equations.
    - If the vectorfield does not include variational equations (or if you
      want second variationals) give hnd>0 and equal to the differentiation
      step to be used for the numerical computation of the stm.

    h,hmin,hmax,tol are the usual parameters for the RKF integrator.
    Doesn't matter the sign of the initial step h. Integration always goes
    from ti to tf determining this way the initial (pos or neg) step.
    '''
    ndia=xi.shape
    ndim=max(ndia) #huh???
    if (abs(ti-tf) < 1.e-12):
        xf=xi
        stm=np.eye(N=ndim)
        return
    h=abs(h)
    if hnd==0:   #Computation of STM using variational equations #TODO not tested yet
        xa=np.zeros(shape=[1,ndim*(ndim+1)])
        xa[0:ndim]=xi
        for i in range(ndim+1,ndim*(ndim+1),ndim+1):
            print(i)
            xa[i]=1.e0
        [_,xfa]=propTiTf(ti,xa,tf,vfield,h,hmin,hmax,tol,params)
        xfa = xfa[-1]
        xf=xfa[0:ndim]
        stm=np.zeros((ndim,ndim))
        for i in range(ndim):
            stm[:,i]=xfa[ndim*i+1:ndim*(i+1)] #TODO idem, error probs
    else:  #Numerical comptuation of STM
        [_,xf]=propTiTf(ti,xi,tf,vfield,h,hmin,hmax,tol,params)
        xf = xf[-1]
        stm=np.zeros((ndim,ndim))
        for k in range(ndim):
            xa=copy.deepcopy(xi)
            xa[k]=xi[k]+hnd
            [_,xfp]=propTiTf(ti,xa,tf,vfield,h,hmin,hmax,tol,params)
            xfp = xfp[-1]
            xa=copy.deepcopy(xi)
            xa[k]=xi[k]-hnd
            [_,xfm]=propTiTf(ti,xa,tf,vfield,h,hmin,hmax,tol,params)
            xfm = xfm[-1]
            #print("xfp",xfp)
            #print("xfm",xfm)
            stm[:,k]=(xfp-xfm)/(2*hnd)#TODO saltarà error crec
    return [xf, stm]

def propTallSec(ti,xi,h,hmin,hmax,tol,vfield,gdgsec,ntall,params):
    '''
    Propagates the inital state xi at time ti with stopping condition g(x)=0
    after ntall cuts. Returns final state or trajectory.

    NOTE: usually ntall > 0, but also one can give ntall=0 and the Newton
          procedure starts from ti, xi to find tf, xf with g(xf)=0.

    Function g(x) must be provided with its row gradient vector in gdgsec.
    This is, [g,dg]=gdsec(x).

    h,hmin,hmax,tol,vfield: usual rk45f input parameters. In particular the
    sign of h determines the sense of propagation and vfield the vectorfield.

    output tf and xf are final time and state with g(xf)=0
'''
    t=ti
    x=xi
    if ntall > 0:
        [t,x,h,_]=rk45f(t,x,h,hmin,hmax,tol,vfield,params)
        [gv,_]=gdgsec(x)
    #-----------------------------------------------------------
    if ntall>0:
        tf = []
        xf = []
        for nt in range(ntall):
            gva=gv
            while gv*gva>0:
                [t,x,h,_]=rk45f(t,x,h,hmin,hmax,tol,vfield,params)
                gva=gv
                [gv,dgv]=gdgsec(x)
    #------------------------------------------------------
    if ntall == 0:
        [gv,dgv]=gdgsec(x)
    ndimv=len(dgv)
    while abs(gv) > 1.e-12:
        df=vfield(x,params)
        h=-gv/np.dot(dgv,df[0:ndimv])
        [t,x,_,_]=rk45f(t,x,h,hmin,hmax,tol,vfield,params)
        [gv,dgv]=gdgsec(x)
    tf=t
    xf=x
    xf = np.array(xf)
    tf = np.array(tf)
    #print(xf.shape,"\n",xf)
    return [tf, xf]


def calDifPoincNum(ti,xi,nd,ntall,gdgSec,pasd,h,hmin,hmax,tol,vfield,params):
    '''
    Numerical computation of the differential of a Poincare map, both in coord
    and time. This is, the Poincare map is [yf tf]=P(xi,ti), where tf and
    yf are the landing time and state in a section secgdg g(x)=0 after ntall
    cuts.
    At output dPoXN contains the differential of yf wrt xi and dPoXt the
    differential of tf wrt xi. (i.e. differential wrt ti which could be
    interesting for non autonomous vfields is not given).

    nd is the dimension of the differential matrix, which will be computed
    considering variations of the first nd components of xi with a numerical
    differentiation step pasd.

    Function g(x) must be provided with its row gradient vector in gdgsec.
    This is, [g,dg]=gdsec(x).

    h,hmin,hmax,tol,vfield: usual rk45f input parameters. In particular the
    sign of h determines the sense of propagation and vfield the vectorfield.
    '''
    dPoXN=np.zeros((nd,nd))
    dPoTN=np.zeros((nd))
    for i in range(nd):
        xa=copy.deepcopy(xi)
        xa[i]+=pasd
        [tfp, xfap] = propTallSec(ti,xa,h,hmin,hmax,tol,vfield,gdgSec,ntall,params)
        xa=copy.deepcopy(xi)
        xa[i]-=pasd
        [tfm, xfam] = propTallSec(ti,xa,h,hmin,hmax,tol,vfield,gdgSec,ntall,params)
        #print('x,y,z    i=',i,' final plus: ',xfap[:3])# Some prints for a common case
        #print('xd,yd,zd i=',i,' final plus: ',xfap[3:6])
        #print('x,y,z    i=',i,' final minus: ',xfam[:3])
        #print('xd,yd,zd i=',i,' final minus: ',xfam[3:6])
        xfap_sec = xfap[0:nd]
        xfam_sec = xfam[0:nd]
        dPoXN[:,i]=(xfap_sec-xfam_sec)/(2.e0*pasd)
        dPoTN[i]=(tfp-tfm)/(2.e0*pasd)
    return [dPoXN, dPoTN]

def fdfArcToSection(ti,xi,h,hmin,hmax,tol,vfield,gdgSec,ntall,hnd,params):
    '''
    It computes the Poincare map and its differential from the state (ti,xi)
    integrated in vectorfiled vfield till the cut number ntall to the section
    given in gdgSec. Output (tf, xf) is the state on gdgSec and dPoX, dPoT the
    differential of the Poincare map [tf, xf]=Po(ti,xi) wrt xi.

    The computation of dPo's can be either obtained from variational
    equations (when these equations are also in the vectorfield, give hnd=0)
    or by numerical differentiation, if the vectorfiled doesn't contain them
    (in this later case give hnd>0 and equal to the differentiation step).

    h,hmin,hmax,tol are the usual parameters for the RKF integrator.
    '''
    ndia=xi.shape
    ndim=max(ndia)
    if hnd==0:   #Computation of dPo using variational equations
        xia=np.zeros((1,ndim*(ndim+1)))
        xia[1:ndim]=xi
        for i in range(ndim+1,ndim*(ndim+1),ndim+1):
            xia[i]=1.e0
        [tf, xfa] = propTallSec(ti,xia,h,hmin,hmax,tol,vfield,gdgSec,ntall,params)
        stm=np.zeros((ndim,ndim))
        for i in range(ndim):
            stm[:,i]=xfa[ndim*i+1:ndim*(i+1)]
        [dPoX,dPoT]= calDifPoinc(tf,xfa,stm,gdgSec,vfield,params) #TODO WARN not yet implemented
        xf=xfa[1:ndim]
    else:  #Numerical comptuation of dPo
        [tf, xf] = propTallSec(ti,xi,h,hmin,hmax,tol,vfield,gdgSec,ntall,params)
        [dPoX,dPoT] = calDifPoincNum(ti,xi,ndim,ntall,gdgSec,hnd,h,hmin,hmax,tol,vfield,params)
    return [tf,xf,dPoX,dPoT]

def cal_FDF_for_PO(tv,xv,pt,ntall,h,hmin,hmax,tol,vfield,gdgSec,hnd,params):
    '''
    Computes function F (where F=0 means periodic orbit) and its differential
    DF (in fact it is the diff of -F), such that solving DF*DX=F we obtain
    the correction step in an iterative process of computing a periodic orbit.

    tv(narcs) (row vector!) and xv(narcs,ndim) is the table of times and
    respective states where we evaluate F. If narcs>1, a parallel shooting is
    performed from t(i) to t(i+1) for i=1...narcs-1. When narcs=1, there are
    no parallel shooting arcs.
    The final arc (the only one when narcs=1) is a Poincare map that proceeds
    from tv(i) xv(i,:) with i=narcs, till the section gdgSec, stopping after
    ntall cuts. The epoch of this final cut and the corresponding state is
    returned in tf,xf.

    ** Description of F:
     The function F=F(x) where x=[xv(1,:)'xv(2,:)',...,xv(k,:)'] and so the
     dimension of x is k*ndim. F is a column vector where its (k-1)*ndim
     first rows represent the matching errors of the parallel shooting
     (the first k-1 blocs of dimension ndim). It is followed by a block of
     rows of dimension ndim giving the matching error of the Poincare map of
     xv(k,:) with xv(1,:) (i.e. the closing condition xv(1,:)=P(xv(k,:) at
     the ntall encounter with gdgSec). This is F when pt=0.
     When pt>0, F includes an aditional row (row number k*ndim+1) which
     contains pt-tprop where tprop is the total time of the propagation
     from tv(1) till the epoch when the last arc ends in the section. This
     is, if pt represents the desired period of the orbit, F(k*ndim+1)
     measures the deviation of the propagation time wrt this period.
     Note that, as tv(*) are assumed fixed, then F(k*ndim+1) is a
     function of the last state xv(k,:).

     DF is de Jacobian of -F with respect to x (defined above). It has
     dimension (k*ndim,k*ndim) when pt=0 and (k*ndim+1,k*ndim) when pt>0.

     tf,xf is the final stopping time and state on the section gdgSec.

    h,hmin,hmax,tol are the usual parameters for the RKF integrator of the
    vectorfield vfield.

    hnd is used to compute stm and dPo from variational equations (hnd=0)
    or numericall (hnd=derivatve_step). See comments inside the functions
    fdfArcFromTiToTf or fdfArcToSection.
    '''
    #print("cal_FDF_for_PO")
    narcs,ndim=xv.shape
    naux=narcs*ndim
    if pt > 0:
        naux=naux+1
    F=np.zeros(shape=(naux))
    DF=np.zeros(shape=(naux,narcs*ndim))
    if narcs>0:
        for k in range(narcs-1): #fa narcs-1 iteracions
            ti=tv[k]
            xi=xv[k]
            tf=tv[k+1]
            [xf, DFk] = fdfArcFromTiToTf(ti,xi,tf,h,hmin,hmax,tol,vfield,hnd,params)
            llr_from = k*ndim
            llr_to = (k+1)*ndim
            llc_from = k*ndim
            llc_to = (k+2)*ndim
            F[llr_from:llr_to]=xv[k+1]-xf
            #DF és una block matrix 2x2. DF[0][0] és DFk, i DF[0][1] és -Id
            DF[llr_from:llr_to,llc_from:llc_to]=np.hstack((DFk,-np.eye(N=ndim)))
    ti=tv[narcs-1]
    xi=xv[narcs-1]
    [tf,xf,dPoX,dPoT] = fdfArcToSection(ti,xi,h,hmin,hmax,tol,vfield,gdgSec,ntall,hnd,params)
    lloc=(narcs-1)*ndim
    llr_from=lloc
    llr_to=lloc+ndim
    F[llr_from:llr_to]=xv[0]-xf
    #DF és una block matrix 2x2. DF[1][0] és -Id, i DF[1][1] és dPoX
    DF[llr_from:llr_to,llr_from:llr_to]=dPoX
    DF[llr_from:llr_to,0:ndim]=DF[llr_from:llr_to,0:ndim]-np.eye(N=ndim)
    if pt > 0: #TODO pt probably not working yet
        F[narcs*ndim]=pt-tf+tv[0]
        DF[narcs*ndim,llr]=dPoT
    return [F,DF,tf,xf]

def rtbp(x,params):
    mu=0.1
    r12= (x[0]-mu+1)**2 + x[1]**2 + x[2]**2 # r12: square of distance to P
    r22= (x[0]-mu)**2 + x[1]**2 + x[2]**2   # r22: square of distance to S
    r13=r12*np.sqrt(r12)
    r23=r22*np.sqrt(r22)

    Ox = x[0] - ((1-mu)*(x[0]-mu)/r23 + mu*(x[0]-mu+1)/r13)
    Oy = x[1] - ((1-mu)* x[1]    /r23 + mu* x[1]      /r13)
    Oz =      - ((1-mu)* x[2]    /r23 + mu* x[2]      /r13)

    xdot = np.zeros(6)
    xdot[0] = x[3]
    xdot[1] = x[4]
    xdot[2] = x[5]

    xdot[3] = 2*x[4] + Ox
    xdot[4] =-2*x[3] + Oy
    xdot[5] = Oz

    return xdot

def tableToPlotPO(tv,xv,tf,nvvf,h,hmin,hmax,tol,vfield,params):
    '''
    Fa taula per fer plot d'orbita i veure com evoluciona el procediment
    de correccio o com es el resultat final del calcul OP.
    tv, xv t i estats (tot en files). tf temps de punt final (usualment sobre
    la seccio final)
    nvvf nobre de variables que te vfield (per si s'usa camp amb variacionals
    que nomes agafi equacions d'estat pos+vel)
    TODO: * Posar flag per decidir si matching points s'inclouen dues vegades
            (per exemple per veure discont en el proces de correccio OP) o be
            guardem a partir de la segona component de sortida de propTiTf.
          * Tambe es podria fer que la h no comences sempre del mateix valor
            a cada arc (normalment potser petit), agafant el valor que tenia
            al final de l'arc anterior e.d usar h=tt(end)-tt(end-1) per comencar
            la propagacio del nou arc.
          * Potser millor en lloc d'usar tf, incloure el darrer arc (el que
            esta fora del for) amb propagacio de Poincare fins seccio, es
            a dir usant propTallSec en lloc de propTiTf
    '''
    [narcs,ndim]=xv.shape
    xi=np.zeros(nvvf)
    tabTX=[]
    for k in range(narcs-1):
        ti=tv[k];
        xi[0:ndim]=xv[k]
        tfa=tv[k+1];
        [tt, xt] = propTiTf(ti,xi,tfa,vfield,h,hmin,hmax,tol,params)
        tabTX.append([tt,xt[:,:ndim]])
    ti=tv[narcs-1]
    xi[0:ndim]=xv[narcs-1]
    [tt, xt] = propTiTf(ti,xi,tf,vfield,h,hmin,hmax,tol,params);
    tabTX.append([tt,xt[:,:ndim]])
    return tabTX

def compute_op(params: list, xi: list):
    '''
    Input:
        params: galparams list
        xi: xinicial
    '''
    print("\ncompute_op")
    np.set_printoptions(precision=5)
    np.set_printoptions(linewidth=np.inf)
    hmin=1.e-4
    hmax=1.e0
    tol=1.e-10
    h=0.01
    ti=0
    gdgSec = lambda xvec : [xvec[1],[0,1,0,0,0,0]] #g(xvec)=y, dg(xvec)=[0,1,0,0,0,0]

    deltaxi = np.array([1e-5,0,-1e-4,0,0,0]) #TODO perquè aquest valor en concret?

    xi += deltaxi

    #fixa les condicions de càlcul
    tolSVD=1.e-5 #mes petit que aixo considero valor singular massa petit

    ntp=1    #nombre d'arcs de tir paral.lel i pas
    htp=1    #delta temps pel tir paralel
    ntall=2  #vigila el nombre de talls fins a seccio que queden !!

    cjv=0
    pt=0
    #pt=2.45 Tria periode o deixa els 0 si vols iteracio norma minima.
    #%cji=cjrtbp(xi,mu)
    #cjv=3.631 Descomentant aquesta fara target cj, pero pt ha de ser >0 pel lloc a DF !

    hnd=1e-5  #Posa 0 per indicar que vfield te variacionals, 1.e-5 per derivacio num de stm
    if hnd==0: #TODO no implementat encara!
        vfield=CAMPSB_var
    else:
        vfield=CAMPSB

    # Genera una semilla inicial con ntp arcos de parall shoot
    tv=[ti] #llista de temps
    xv=[xi] #llista de posicions
    for i in range(ntp):
        ta=ti+(i+1)*htp
        [tfa, xfa] = propTiTf(ti,xi,ta,vfield,h,hmin,hmax,tol,params) #%0 del finalno tabla
        tfa = tfa[-1]
        xfa = xfa[-1]
        tv.append(tfa) #CONCAT last point
        xv.append(xfa) #CONCAT last point

    print("xv,",xv)
    print("tv,",tv)
    xv = np.array(xv)
    tv = np.array(tv)

    narcs,ndim=xv.shape #num d'arcs d'orb que s'integren i dim del camp vect
    print("narcs,",narcs,"ndim,",ndim)
    print("\n")

    for i in range(10): #TODO hardcode???
        [F,DF,tf,xf]=cal_FDF_for_PO(tv,xv,pt,ntall,h,hmin,hmax,tol,vfield,gdgSec,hnd,params)
        #print(i,"F")
        #print(F)
        #print(i,"DF")
        #print(DF)

        ############################################################
        if cjv > 0:  #canviem equacio de target period per equacio target Ct Jac
            [cjc, dcj]= cjdcj(xv[narcs-1]) #TODO potencialment això no funciona pels index
            llr_from=(narcs-1)*ndim+1
            llr_to=narcs*ndim
            F[narcs*ndim]=cjv-cjc
            DF[narcs*ndim,llr_from:llr_to]=dcj

        #Correccio usant SVD, corr=V*inv(S)*U'*F amb svd's "prou grans"
        [U,diagS,V] = np.linalg.svd(DF)
        S = np.diag(diagS) #keep it consistent with MatLab
        nsf,nsc = DF.shape

        #count nonzero values of diagonal
        nvsg = int(np.sum([1 if x>tolSVD else 0 for x in diagS]))

        U[:,nvsg+1:nsf]=[] #TODO sembla que no passa?
        SA=1./diagS[0:nvsg]
        V[:,nvsg+1:nsc]=[]  #TODO sembla que no passa?
        #print("U\n",U)
        #print("SA\n",SA)
        #print("V\n",V)
        for j in range(nvsg):
            U[:,j]=U[:,j]*SA[j]
        #print("U\n",U)
        UxF = np.matmul(U.transpose(),F) #hauria de ser 1x12
        #print("A",UxF)

        #print("UxF",UxF)
        dxv = np.matmul(V.transpose(), UxF) #hauria de ser 1x12
        #print("dxv",dxv)
        #print("xv correction iter ",i)
        #print("dxv",dxv)
        for k in range(narcs):
            xv[k] += dxv[(k)*(ndim):(k+1)*(ndim)] #TODO revisar això que segur que ho puc deixar millor
        print(xv)

        #%----- dibuix
    tabTX = tableToPlotPO(tv,xv,tf,ndim,h,hmin,hmax,tol,vfield,params)
    print(len(tabTX[1][0]),len(tabTX[1][1][0]))
    temps = np.hstack((tabTX[0][0],tabTX[1][0]))
    pos = np.vstack((tabTX[0][1],tabTX[1][1])).transpose()
    plt.scatter(pos[0],pos[1],c=temps)
    plt.show()
        #hold on
