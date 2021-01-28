using Plots; gr()
using Distributions
using Statistics
using LinearAlgebra
using JSON

# This is a helper function that will be useful going forward...
# It constructs SymmetricGroup( A ), and the elements of A can be anything
# note: S_k = permGroup( collect(1:k) )
function permGroup(A::Union{Set,Array})
    unique!(A)
    P = Array{eltype(A),1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(P, head)
        end
    end
    continuePerm(eltype(A)[], A)
    return P
end


permGroup(collect(1:12))

# this constructs a matrix for a given permutation
# inSnWithn
function PMatrix(τ::Array; inSnWithn=nothing)
    if inSnWithn==nothing
        p = zeros(Int16, length(τ),length(τ))
    else
        p = zeros(Int16, inSnWithn,inSnWithn)
    end
    for i in 1:length(τ)
        p[ i , τ[i] ] = 1
    end
    return p
end

# helper function for checking if it is valid to say X < Y
# note: differs from julia's built in issorted which seems to use partial orders
#   ex: issorted([1,2,2,3]) returns true
function isordered(X::Array{T},Y::Array{T}) where T <: Number
    if all(map(t-> isless(t[1],t[2]), Base.product(X,Y)))
        return true
    else
        return false
    end
end

# this tests if an array is strictly ordered by <, which is defined for various things
function strictOrder(S)
    for i in 1:(length(S)-1)
        if !(S[i] < S[i+1])
            return false
        end
    end
    return true
end

# Assumes we pass
function totVar(A)
    n = size(A)[1]
    v = 0
    for τ in permGroup( collect(1:n) )
        v += abs( reduce(*,[ A[i,τ[i]] for i in 1:n ]) )
    end
    return v
end


# This makes dealing with data quite a bit easier,
#   though the use of symbols can cause difficulties & headaches
isoDict = Dict([  "Australia" => :AUS,
                "Basque Country" => :ESP,
                "Brazil" => :BRA,
                "Fiji" => :FJI,
                "France" => :FRA,
                "Hawaii" => :USA,
                "Indonesia" => :IDN,
                "Italy" => :ITA,
                "Japan" => :JPN,
                "New Zealand" => :NZL,
                "Portugal" => :PRT,
                "South Africa" => :ZAF,
                "Spain" => :ESP,
                "United States" => :USA ])

isoS = sort!(unique(collect(values(isoDict)))) # For convenience

# Loading data
data = Dict()
data = JSON.parse(open("Data/CombinedCountries/CleanAllDataCC.txt", "r"))  # parse and transform data

waves = []
for wid in keys(data)
    if data[wid]["nJudOrigs"] == 5 & data[wid]["nSubScos"] == 5
        origs = unique(data[wid]["subScoOrig"])
        matchIndicator = (data[wid]["athOrig"] in origs)
        labeledScos = Dict([isoDict[origin] => Float16[] for origin in origs])
        origScoPairs = collect(zip(data[wid]["subScoOrig"],data[wid]["subSco"]))

        labeledScosBinary = Dict([:Match => Float16[], :NoMatch => Float16[] ])

        for p in origScoPairs
            # push!( array of judge scores from country p[1], score=p[2] )
            push!(labeledScos[ isoDict[p[1]] ], p[2])
            if p[1] == data[wid]["athOrig"]
                push!(labeledScosBinary[:Match], p[2])
            else
                push!(labeledScosBinary[:NoMatch], p[2])
            end
        end

        x = (   id=wid,
                evtYear=data[wid]["evtYear"],
                evtOrig=isoDict[data[wid]["evtOrig"]],
                evtName=data[wid]["evtName"],
                evtId=data[wid]["evtId"],
                rnd=data[wid]["rnd"],
                rndId=data[wid]["rndId"],
                heat=data[wid]["heat"],
                heatId=data[wid]["heatId"],
                athName=data[wid]["athName"],
                athId=data[wid]["athId"],
                athOrig=isoDict[data[wid]["athOrig"]],
                currentPoints=data[wid]["currentPoints"],
                endingPoints=data[wid]["endingPoints"],
                panel=labeledScos,
                panelBinary=labeledScosBinary,
                subScos=data[wid]["subSco"],
                subScoOrigs=map(x->isoDict[x], data[wid]["subScoOrig"]),
                panelOrigs=Set(map(x->isoDict[x], data[wid]["subScoOrig"])),
                match=matchIndicator )

        push!(waves, x)
    end
end

# Generally Used, and helpful
hasMatch = filter(w->w.match, waves)


Surfer = BRA

[AUS, BRA, AUS, USA, NZT]

Match < NoMatch
NoMatch < Match


#-----------------------------------------------------
EVTS = sort(unique(map(x->x.evtName, waves)))
TotVarEvts = zeros(length(EVTS))
U = [1/2 1/2; 1/2 1/2]
for (i,vt) in enumerate(EVTS)
    thisEventPanels =  map(x-> x.panelBinary, filter(x->x.evtName==vt, hasMatch))
    D = zeros(2,2)
    noOrder = 0
    for panel in thisEventPanels
        if isordered(panel[:Match],panel[:NoMatch])
            D += [1 0; 0 1]
        elseif isordered(panel[:NoMatch], panel[:Match])
            D += [0 1; 1 0]
        else
            noOrder +=1
            D += U
        end
    end
    D = length(thisEventPanels)^-1 * D
    println(length(thisEventPanels), " panels at ", vt, " and ", noOrder, " have no order")
    TotVarEvts[ i ] = totVar(D-U)
end
bar(EVTS, TotVarEvts)



#---------------------------------------------
HEATS = sort!(unique(map(x->x.heatId, waves)))
TotVarHts = zeros(length(HEATS))
U = [1/2 1/2; 1/2 1/2]
for (i,ht) in enumerate(HEATS)
    thisHeatPanels = map(x-> x.panelBinary, filter(x->x.heatId==ht, hasMatch))
    D = zeros(2,2)
    noOrder = 0
    for panel in thisHeatPanels
        if isordered(panel[:Match],panel[:NoMatch])
            D += [1 0; 0 1]
        elseif isordered(panel[:NoMatch], panel[:Match])
            D += [0 1; 1 0]
        else
            D += U
            noOrder += 1
        end
    end
    D = length(thisHeatPanels)^-1 * D
    if isapprox(D,U)
        TotVarHts[i] = 0
    else
        #println(length(thisHeatPanels), " panels in Heat ", ht, " and ", noOrder, " have no order")
        t = totVar(D-U)
        TotVarHts[i] = t
    end
end
bar(TotVarHts)
histogram(TotVarHts)

63746
63762
68026
72037
72707---
72714---
77674
78716---


#---------------------------------------------
SURFERS = sort(unique(map(x->x.athName, waves)))
surferToIndex = Dict([sfr => i for (i,sfr) in enumerate(SURFERS)])
nS = length(SURFERS)
totalSurfer = zeros(nS,nS)
for ht in HEATS
    thisHeat= filter(x->x.heatId==ht, waves)
    athletes = sort(unique(map(x->x.athName, thisHeat)))
    finalScores = []

    athsNotInHeat = setdiff(1:nS, map(x->surferToIndex[x], athletes))
    totalSurfer[athsNotInHeat, athsNotInHeat] .+= 1 / length(athsNotInHeat)

    for ath in athletes
        scores = map( x-> mean( sort(x.subScos)[2:4] ), filter(x->x.athName==ath, thisHeat))
        if length(scores) >=2
            heatScore = sum(sort(scores)[1:2])
        else
            heatScore = sum(scores)
        end
        push!(finalScores, (heatScore,ath))
    end

    sort!(finalScores)
    for (i,ath) in enumerate(athletes)
        totalSurfer[ surferToIndex[ath] ,surferToIndex[finalScores[i][2]] ] += 1
    end
end
heatmap(totalSurfer)


#------------------------------------------------------
Wbin = []
U = [1/2 1/2; 1/2 1/2]
D = zeros(2,2)
noOrder = 0
for w in hasMatch
    if isordered(w.panelBinary[:Match],w.panelBinary[:NoMatch])
        D += [1 0; 0 1]
        push!(Wbin, [1 0; 0 1])
    elseif isordered(w.panelBinary[:NoMatch], w.panelBinary[:Match])
        D += [0 1; 1 0]
        push!(Wbin, [0 1; 1 0])
    else
        D += U
        push!(Wbin, U)
        noOrder += 1
    end
end


#--------------------------------
U = [1/2 1/2; 1/2 1/2]
for country in isoS
    BinaryPanelByAthOrig = map(w-> w.panelBinary, filter(x-> x.athOrig == country, hasMatch))
    noOrderBinary = 0
    yesOrderBinary = 0
    DBinary = zeros(2,2)
    for pb in BinaryPanelByAthOrig
        if isordered(pb[:Match],pb[:NoMatch])
            DBinary += [1 0; 0 1]
            yesOrderBinary += 1
        elseif isordered(pb[:NoMatch], pb[:Match])
            DBinary += [0 1; 1 0]
            yesOrderBinary += 1
        else
            DBinary += U
            noOrderBinary += 1
        end
    end
    println(country)
    println(DBinary)
    DBinary *= length(BinaryPanelByAthOrig)^-1
    println("$yesOrderBinary waves were consistent with a total order {:Match, :NoMatch}")
    println("$noOrderBinary waves were not consistent with any total order on {:Match, :NoMatch}")
    if isapprox(DBinary-U, zeros(2,2))
        println("Total Variation to uniform is 0")
    else
        println(totVar(DBinary-U))
    end
    println("----------------------")
end


#-----------------------------------------
allPanelCompositions = unique( map(x->x.panelOrigs, waves))
judOrigs = sort(collect(∪(allPanelCompositions...)))
permGroup(judOrigs)
U = 1/7 * ones(7,7)
NatToIndex = Dict([orig => i for (i,orig) in enumerate(judOrigs) ])
W = []
consWave = []
totalConsistentPanels = 0
for panelComp in allPanelCompositions
    consistentPanels = 0
    notOnPanel = map(x->NatToIndex[x], setdiff(judOrigs,panelComp) )
    unifM = zeros(7,7)
    unifM[notOnPanel,notOnPanel] .= 1 / length(notOnPanel)
    e = sort(collect(panelComp))
    for w in filter(x -> x.panelOrigs == panelComp, waves)
        waveMatrix = zeros(7,7)
        for ord in permGroup(collect(panelComp))
            S = Iterators.product([w.panel[c] for c in ord]...)
            if all(map(strictOrder, S) )
                for (i,ei) in enumerate(e)
                    waveMatrix[ NatToIndex[ei], NatToIndex[ ord[i]] ] = 1
                end
                waveMatrix += unifM
                consistentPanels +=1
                push!(consWave, w)
                push!(W, waveMatrix)
            end
        end
    end
    totalConsistentPanels += consistentPanels
    println(panelComp)
    println(consistentPanels)
end
totVar(totalConsistentPanels^-1 * sum(W) - U)
