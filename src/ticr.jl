"""
    ticr!(net, D::DataFrame, optimizeBL::Bool)
    ticr!(net, D::DataCF,    optimizeBL::Bool)
    ticr(D::DataCF)

Goodness-of-fit test for the adequacy of the multispecies network coalescent,
to see if a given population or species network explains the
quartet concordance factor data adequately (see Stenz et al 2015
and [addendum](http://www.stat.wisc.edu/~ane/publis/2015Stenz_TICR_addendum.pdf)
for the method on trees)
Obviously the acronym TICR (Tree Incongruence Checking with R) becomes outdated:
with a method extended to networks and implemented in Julia. Oh well :smiley:.

The tree / network needs to have branch lengths in coalescent units,
must be fully resolved, and must be of level 1.

The model assumes a Dirichlet distribution for the observed quartet concordance
factor, with concentration parameter estimated from the data. An outlier p-value
is calculated for each four-taxon set. Four-taxon sets are then binned into
categories according to their p-values: `0-0.01`, `0.01-0.05`, `0.05-0.10`, and `0.10-1`.
Finally, a one-sided goodness-of-fit test is performed on these binned
frequencies to determine if they depart from an
expected 5% p-values being below 0.05 (versus more than 5%).
Note: this departs from Stenz et al. (2015), who tested for any kind of departure
from the expected frequencies (0.01, 0.04, 0.05, 0.90) across all 4 bins.

- The first version takes a DataFrame object where each row corresponds to a given
  four-taxon set. The DataFrame is modified by having an additional another column
  containing the p-values corresponding to each four-taxon set.
- The second version takes a DataCF object and modifies it by updating
  the expected concordance factors stored in that object.
- The last version (which all others call) assumes that the expected concordance
  factors in the DataCF object are correctly calculated from the test network.

`optimizeBL`: when false, the loglik field of `net` is updated;
when `true`, a copy of `net` with updated branch lengths (in coalescent units)
and update loglik is returned.

`minimum`: when true, the outlier p-value per four-taxon set corresponds to the minimum
p-value of the three quartet concordance factors (as opposed to only using the p-value
for the major CF). Default is false.

output:

- p-value of the overall goodness-of-fit test
- test statistic (z value)
- binned frequencies of outlier p-values
- value of the concentration parameter α
- value of the pseudo likelihood (optimized at α)
- a vector of outlier p-values, one for each four-taxon set
- network (first and second versions):
  `net` with loglik field updated if `optimizeBL` is false;
   copy of `net` with optimized branch lengths and loglik if `optimizeBL` is true

# References

NWM Stenz, B Larget, DA Baum and C Ané (2015).
Exploring tree-like and non-tree-like patterns using genome sequences:
An example using the inbreeding plant species *Arabidopsis thaliana* (L.) Heynh.
Systematic Biology, 64(5):809-823. doi: 10.1093/sysbio/syv039
"""
function ticr!(net, D::DataFrame, optimizeBL::Bool; minimum=false::Bool)
    d = readTableCF(D);
    res = ticr!(net, d, optimizeBL,minimum); # order of value in results "res":
    D[:p_value] = res[6]             # overallpval, teststat, counts, alpha, pseudolik, pval, net
    return res
end

function ticr!(net, D::DataCF, optimizeBL::Bool,minimum::Bool)
    if optimizeBL
        net = topologyMaxQPseudolik!(net,D);
    else
        topologyQPseudolik!(net,D);
    end
    if minimum
        res = ticr_min(D);
    else
        res = ticr(D);
    end
    return (res..., net) # (overallpval, teststat, counts, alpha, pseudolik, pval, net)
end

@doc (@doc ticr!) ticr
function ticr(D::DataCF)
    res_alpha = ticr_optimalpha(D)
    alpha = res_alpha[2][1]
    pseudolik = res_alpha[1]
    N = D.numQuartets
    pval = Float64[]
    for i in 1:N
        phat = D.quartet[i].obsCF
        p = D.quartet[i].qnet.expCF
        max = findmax(p)
        p_max = max[1]
        max_idx = max[2]
        p_max_hat = getindex(phat,max_idx)
        p_sort = sort(D.quartet[i].qnet.expCF)
        abs(p_max-p_sort[end-1]) > 1e-6 || warn("Check the network for major quartet")
        d = abs(p_max_hat - p_max)
        temp = [1-(1-p_max)*alpha/2, 0.0]
        shapeAdd = maximum(temp)
        ipval = StatsFuns.betacdf(alpha*p_max+shapeAdd, alpha*(1-p_max)+2*shapeAdd, p_max-d)+
        StatsFuns.betaccdf(alpha*p_max+shapeAdd, alpha*(1-p_max)+2*shapeAdd, p_max+d)
        push!(pval,ipval)
    end
    pcat = CategoricalArrays.cut(pval,[0, 0.01, 0.05, 0.1, 1]) # CategoricalArrays is required by DataFrames
    count = StatsBase.countmap(pcat)
    c = Float64[]
    e = [0.01,0.04,0.05,0.90]*N
    interval = ["[0.0, 0.01)","[0.01, 0.05)","[0.05, 0.1)","[0.1, 1.0)"]
    for i in interval
        if haskey(count, i)
            push!(c,count[i])
        else
            push!(c,0.0)
        end
    end
    # in Stenz et al (2015) and in phylolm implementation:
    # teststat = 0.0 # chi-square test statistics
    # for i in 1:length(c) teststat += (c[i]-e[i])^2 / e[i]; end
    # overallpval = StatsFuns.chisqccdf(3,teststat)
    # instead: reduce to [0,.05) and [.05,1.0] and do a one-sided test
    teststat = ((c[1]+c[2])/N - 0.05)/sqrt(0.0475/N) # z-value. 0.0475 = 0.05 * (1-0.05)
    overallpval = StatsFuns.normccdf(teststat) # one-sided: P(Z > z)
    return (overallpval, teststat, count, alpha, pseudolik, pval)
end

"""
    ticr_optimalpha(D::DataCF)

Find the concentration parameter α by maximizing the pseudo-log-likelihood
of observed quartet concordance factors. The model assumes a Dirichlet distribution
with mean equal to the expected concordance factors calculated from a
phylogenetic network (under ILS and reticulation). These expected CFs
are assumed to be already calculated, and stored in D.

When calculating the pseudo-log-likelihood, this function will check the
observed concordance factors for any values equal to zero: they cause a problem
because the Dirichlet density is 0 at 0 (for concentrations > 1).
Those 0.0 observed CF values are re-set to the minimum of:
- the minimum of all expected concordance factors, and
- the minimum of all nonzero observed concordance factors.

output:

- maximized pseudo-loglikelihood
- value of α where the pseudo-loglikelihood is maximized
- return code of the optimization

The optimization uses NLOpt, with the `:LN_BOBYQA` method.
Optional arguments can tune the optimization differently:
`NLoptMethod`, `xtol_rel` (1e-6 by default),
starting α value `x_start` (1.0 by default).
"""
function ticr_optimalpha(D::DataCF; x_start=1.0::Float64,
        NLoptMethod=:LN_BOBYQA::Symbol, xtol_rel=1e-6::Float64)
    M = D.numQuartets
    logCFtilde = 0.0
    minExpCF = minimum([minimum(q.qnet.expCF) for q in D.quartet])
    minNonZeroObsCF = minimum([minimum(filter(!iszero,q.obsCF)) for q in D.quartet])
    minObsCF = min(minExpCF,minNonZeroObsCF)
    for q in D.quartet
        q.obsCF[find(q.obsCF .== 0)] = minObsCF
    end
    for i in 1:M
        lobsCF = log.(D.quartet[i].obsCF)
        p = D.quartet[i].qnet.expCF
        for i in 1:length(lobsCF)
            logCFtilde += lobsCF[i]*p[i]
        end
    end
    function obj(x::Vector, grad::Vector)
        a = x[1]
        sum = 0
        logCFbar = 0.0
        for i in 1:M
            lobsCF = log.(D.quartet[i].obsCF)
            p = D.quartet[i].qnet.expCF
            p_min = minimum(p)
            if a >= 1/p_min
                b = 0
            else
                b = 1 - a * p_min
            end
            logCFbar += mean(lobsCF)*(1-b)
            sum += lgamma(a+3*b) - lgamma(a*p[1]+b) - lgamma(a*p[2]+b) - lgamma(a*p[3]+b)
        end
        PseudologLik = sum + a*logCFtilde-3*logCFbar
        return PseudologLik
    end
    opt = NLopt.Opt(NLoptMethod,1)
    NLopt.lower_bounds!(opt, 0)
    NLopt.upper_bounds!(opt, 10^5)
    NLopt.maxeval!(opt,1000)
    NLopt.xtol_rel!(opt,xtol_rel)
    NLopt.max_objective!(opt,obj)
    fmax, xmax, ret = NLopt.optimize(opt,[x_start])
    return fmax, xmax, ret
end


## same as ticr, but instead of choosing the maximum, we do for every CF,
## and take the minimum pvalue
function ticr_min(D::DataCF)
    res_alpha = ticr_optimalpha(D)
    alpha = res_alpha[2][1]
    pseudolik = res_alpha[1]
    N = D.numQuartets
    pval = fill(0.0,N,3)
    for i in 1:N
        phat = D.quartet[i].obsCF
        p = D.quartet[i].qnet.expCF
        for j in 1:3 ##we now do ticr for every CF
            p_max = p[j]
            p_max_hat = phat[j]
            p_sort = sort(D.quartet[i].qnet.expCF)
            abs(p_max-p_sort[end-1]) > 1e-6 || warn("Check the network for major quartet")
            d = abs(p_max_hat - p_max)
            temp = [1-(1-p_max)*alpha/2, 0.0]
            shapeAdd = maximum(temp)
            ipval = StatsFuns.betacdf(alpha*p_max+shapeAdd, alpha*(1-p_max)+2*shapeAdd, p_max-d)+
            StatsFuns.betaccdf(alpha*p_max+shapeAdd, alpha*(1-p_max)+2*shapeAdd, p_max+d)
            pval[i,j] = ipval
        end
    end
    ## for every quartet, we calculate min pvalue
    minpval = fill(0.0,N)
    for i in 1:N
        minpval[i] = minimum(pval[i,:])
    end
    pcat = CategoricalArrays.cut(minpval,[0, 0.01, 0.05, 0.1, 1]) # CategoricalArrays is required by DataFrames
    count = StatsBase.countmap(pcat)
    c = Float64[]
    e = [0.01,0.04,0.05,0.90]*N
    interval = ["[0.0, 0.01)","[0.01, 0.05)","[0.05, 0.1)","[0.1, 1.0)"]
    for i in interval
        if haskey(count, i)
            push!(c,count[i])
        else
            push!(c,0.0)
        end
    end
    # in Stenz et al (2015) and in phylolm implementation:
    # teststat = 0.0 # chi-square test statistics
    # for i in 1:length(c) teststat += (c[i]-e[i])^2 / e[i]; end
    # overallpval = StatsFuns.chisqccdf(3,teststat)
    # instead: reduce to [0,.05) and [.05,1.0] and do a one-sided test
    teststat = ((c[1]+c[2])/N - 0.05)/sqrt(0.0475/N) # z-value. 0.0475 = 0.05 * (1-0.05)
    overallpval = StatsFuns.normccdf(teststat) # one-sided: P(Z > z)

    return (overallpval, teststat, count, alpha, pseudolik, minpval)
end
