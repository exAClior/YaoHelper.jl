"""
This file contains helper functions for creating commonly used states.
"""

export one_state, plus_state, minus_state

function one_state(n)
    return zero_state(n) |> repeat(n,X)
end

function plus_state(n)
   return zero_state(n) |> repeat(n,H)
end

function minus_state(n)
    return zero_state(n) |> repeat(n,H) |> repeat(n,Z)
end
