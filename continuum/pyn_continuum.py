def LEGFIT(spec, minord, maxord):
    """
    Fit Legendre polynomial continua
    """

    def _ret():  return spec

    x = spec['velocity'].copy()
    y = spec['flux'].copy()


    nflag = 0
    nord = maxord
    #
    #Array subscript length and vector.
    #
    nx = x.size
    ix = np.arange(nx, dtype=int32)

    #Form legendre polynomial.
    p = np.zeros([maxord + 1, nx], "float32")
    p[ix] = 1.
    p[ix + nx] = x
    for j in np.arange(2., (maxord)+(1)):
        p(ix + j * nx) = ((2. * j - 1.) * x * p[j - 1,:] - (j - 1) * p[j - 2,:]) / j
    #
    #Begin loop to do fit.
    #
    for nord in np.arange(minord, (maxord)+(1)):

    ## LOOP:
        ncoeff = nord + 1
        #
        #Form alpha and beta matrices.
        #
        beta = np.zeros([nord + 1], "float32")
        alpha = np.zeros([nord + 1, nord + 1], "float32")
        for k in np.arange(0, (nord)+(1)):
            beta(k) = TOTAL(y * p[k,:])
        for k in np.arange(0, (nord)+(1)):
            for j in np.arange(0, (nord)+(1)):
                alpha(j, k) = TOTAL(p[j,:] * p[k,:])
        #
        #Invert alpha matrix ==> error matrix eps.
        #
        eps = INVERT(alpha)
        #
        #Calculate coefficients and fit.
        #
        a = transpose(matrixmultiply(transpose(beta), transpose(eps)))
        yfit = np.zeros([nx], "float32")
        for j in np.arange(0, (nord)+(1)):
            yfit = yfit + a(j) * p[j,:]
        #
        #Calculate chi squared of fit - uniform weighting=1.
        #
        #		sigma = SQRT(TOTAL((y-yfit)^2)/(nx-ncoeff-1))
        chisq = TOTAL((y - yfit) ** 2)
        chi2 = chisq / (nx - ncoeff - 1)
        #
        #Check chi2 against previous chi2 to see if additional term should
        #be kept.  Check probability for "95% confidence".
        #
        # IF nflag EQ 1 THEN  GOTO,OUT
        if nord > minord:
            f = (chisq1 - chisq) / chi2
            fcutoff = FTEST(nx - ncoeff - 1, 0.05)
            if f < fcutoff:
                nord = nord - 1
                nflag = 1
                ## GOTO,;; LOOP  ;back up to top to do it over again
        chisq1 = chisq
    # OUT:
    # 	maxord = maxord < nord



    #
    #Calculate errors on coefficients - not used here since uniform weighting of
    #data, but could be used someday.
    #
    #	siga = FLTARR(maxord+1)
    #	FOR j=0,maxord DO siga(j) = SQRT(chi2*eps(j,j))

    # RETURN

    return _ret()
