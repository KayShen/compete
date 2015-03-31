steepness <- function(M,n_tries = 2000,simulation = FALSE){
    n = nrow(M)
    ds_p = DS(M)
    n_ds_p = NormDS(ds_p)
    ds_d = Dyadic(M)
    n_ds_d = NormDS(ds_d)
    r = n_ds_p[order(n_ds_p)][n:1]
    r2 = n_ds_d[order(n_ds_d)][n:1]
    fit = lm(r~rep(1:n))
    fit2 = lm(r2~rep(1:n))
    slope = coefficients(fit)[2]
    intercept = coefficients(fit)[1]
    slope2 = coefficients(fit2)[2]
    intercept2 = coefficients(fit2)[1]
    steepness_p = -slope
    steepness_d = -slope2
    if (simulation == FALSE){
        plot(r,type = "l", col = "red", ylab = "NormDS", xlab = "Rank order")
        points(r2, type = "l", col = "blue")
        points(r,pch = 15, col = "red")
        points(r2,pch = 16, col = "blue")
        abline(coef = coef(fit), col = "red", lty = 4)
        abline(coef = coef(fit2), col = "blue", lty = 4)
        legend(0.6*n,r[1],c("NormDS based on P","NormDS based on D"),pch = c("l","l"),col = c("red", "blue"))
        sim = simulation(M, n_tries, steepness_p, steepness_d)
        p_p = sim$p_p
        p_d = sim$p_d
        result = list(steepness_p = as.numeric(steepness_p), steepness_d = as.numeric(steepness_d), p_p = p_p, p_d = p_d)
    }
    else{
        result = list(steepness_p = steepness_p, steepness_d = steepness_d)
    }
    return(result)
}

proportion <- function(M){
    M_temp = M + t(M)
    M_temp[M_temp == 0] = 0.01
    M_p = M/M_temp
    return(M_p)
}

DS <- function(M){
    M_p = proportion(M)
    w = rowSums(M_p)
    l = colSums(M_p)
    M_p_weighted = t(w*t(M_p))
    M_p_weighted2 = l*M_p
    w2 = rowSums(M_p_weighted)
    l2 = colSums(M_p_weighted2)
    result = w + w2 - l - l2
    return(result)
}

NormDS <- function(DS){
    n = length(DS)
    result = (DS + n*(n-1)/2)/n
    return(result)
}

Dyadic <- function(M){
    M_p = proportion(M)
    pp = M + t(M)
    pp2 = 1/(pp)
    pp = 1/(1+pp)
    pp[pp2==Inf] = 0
    D = M_p - ((M_p-0.5)*pp)
    M_p = D
    w = rowSums(M_p)
    l = colSums(M_p)
    M_p_weighted = t(w*t(M_p))
    M_p_weighted2 = l*M_p
    w2 = rowSums(M_p_weighted)
    l2 = colSums(M_p_weighted2)
    result = w + w2 - l - l2
    return(result)
}

simulation <- function(M, nTries, slope_p, slope_d){
    Mk = M + t(M)
    total = Mk[upper.tri(Mk)]
    temp = total
    temp2 = total
    steep = rep(NA, nTries)
    steep2 = rep(NA, nTries)
    for (i in 1:nTries){
        for(j in 1:length(total)){
            temp[j] = sample(total[j],1)
        }
        temp2 = total - temp
        Mk[upper.tri(Mk)] = temp
        Mk[lower.tri(Mk)] = temp2
        st = steepness(Mk, 1, TRUE)
        steep[i] = st$steepness_p
        steep2[i] = st$steepness_d
    }
    p_p = sum(steep > slope_p)/nTries
    p_d = sum(steep > slope_d)/nTries
    result = list(p_p = p_p, p_d = p_d)
    return(result)
}
