using EntanglementDist
using SQLite
using Convex

# Since the solver crashes from time to time, one can
# run `./keepGenerating.sh` to restart this script
# automatically after crashing.

### Settings ###
# Database where output is saved
db = SQLite.DB("test.sqlite")
ps = collect(0.4:0.1:1)
δs = collect(0.1:0.05:1)
ϵ = 1e-4 # relative accuracy

# This is the name under which the results are saved
# N.B.
# The program type should be included in the state name as otherwise
# there is no way to tell them apart. `./plotData.jl` assumes
# that the state name is of the following form:
# `{statename}` for PPT program
# `1ext {statename}` for 1 extension
stateName = "pptRelax: Double Werner"
n_copies = 2
#Local dimensions of individual copies:
nA = 2
nB = 2
# local dimension of the desired maximally entangled state
k=2

"""
Generate one copy of the state. Change this to perform numerics on another state!
"""
function generateState(parameter)
  return wernerState(parameter)
  #return rStateCorrPhase(parameter)
end

# No need to change this
entries = length(ps)*length(δs)
times = []
firstRun = true
maxRetries = 20
retries = 0
foundNaN = false

# """
# Generate n copies of the state
# """
# function generateState(parameter, n_copies::Int)
#     return copies(generateState(parameter),  n_copies)
# end

"""
Save found (F, p_succ) to the database
"""
function saveData(F, p_succ, δ, stateName, p, eprF, ϵ)
  q = string("INSERT INTO `pptRelax` (`fidelity`, `p_succ`, `delta`, `state`, `p`, `eprFidelity`, `eps`)
    VALUES (" , F , ", ", p_succ ,", ", δ ,", '", stateName, "', ", p ,", ", eprF, ", ", ϵ ,")")
  SQLite.query(db, q)
end

"""
Check in the database if the combination (state, δ, p_succ) already occurs
"""
function hasBeenComputed(δ, stateName, p, ϵ)
  q = string("SELECT COUNT(*) FROM `pptRelax` WHERE `delta` = ", δ ," AND p = ", p ," AND `state` = '", stateName, "' AND eps <= '", ϵ,"'")
  result = SQLite.query(db, q)
  count = result[1][1].value

  return count > 0
end

"""
Convert seconds to hours, minuts and seconds
"""
function hms(timeleft)
  timeleft = round(Int, timeleft)
  hours = floor(Int, timeleft/3600)
  timeleft = timeleft - hours * 3600

  minutes = floor(Int, timeleft/60)
  timeleft = timeleft - minutes * 60

  return string(hours, "h ", minutes, "m ", timeleft, "s")
end

### Run the programs ###
# sometimes the solved gives an error for unknown reasons
# therefore, retry if an error occurs
while firstRun == true || (foundNaN == true && retries < maxRetries)
  entriesDone = 0
  firstRun = false
  foundNaN = false
  retries = retries + 1

  for p in ps
    for δ in δs
      # start timer
      starttime = time()

      # check if it is already computed
      computed = hasBeenComputed(δ, stateName, p, ϵ)
      
      # only compute data if no entry is available
      if !computed
        state = generateState(p)
        println("computing (δ, p, state) = (", δ, ", ", p, ", ", stateName, ")")

        # Alter the line below to run another program (e.g. PPT or k-ext)
        # The program type should be included in the state name as otherwise
        # there is no way to tell them apart.
        (problem, F, p_succ) = pptRelaxCopies(state, n_copies, nA, nB, k, δ, verbose = false, eps = ϵ, max_iters = 1e7)
        # (problem, F, p_succ) = PPTprogrammeNoTwirling1ExtPermSymOnlySym(state, 2^n_copies, 2^n_copies, n_copies, 2, δ, verbose = true, eps = 1e-4)
        eprF = eprFidelity(generateState(p))
        

        println("found (F, p_succ) = " , (F, p_succ))

        if isnan(F)
          foundNaN = true
        else
          saveData(F, p_succ, δ, stateName, p, real(eprF), ϵ)

          stoptime = time()
          times = [times; stoptime - starttime]
          timeleft = (entries - entriesDone) * mean(times)
          entriesDone = entriesDone + 1
          println("Time to go: ", hms(timeleft), ". computation time: ", hms(stoptime-starttime), " | ", entriesDone, "/" , entries , "\n")
        end
      else
        println("already computed (δ, p, state, ϵ) = (", δ, ", ", p, ", ", stateName, ", ", ϵ ,"). \n")
        entriesDone = entriesDone + 1
      end
    end
  end
end
