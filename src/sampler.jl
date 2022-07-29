function select_event(brates, urates, total_rate)
    threshold = rand()*total_rate
    rate_sum = 0
    overlap_i = 1
    binding = 1
    for i in 1:length(brates)
        rate_sum += brates[i]
        if threshold < rate_sum
            overlap_i = i
            break
        end
        rate_sum += urates[i]
        if threshold < rate_sum
            overlap_i = i
            binding = -1
            break
        end
    end

    return (overlap_i, binding)
end

function sample_trajectory(equation_of_motion, brate_function, urate_function, lambda0,
        Nd0, p, tfinal, dt)
    lambda = lambda0
    lambdas = [lambda]
    Nd = deepcopy(Nd0)
    Ndtotal = sum(Nd)
    Ndtotals = [Ndtotal]
    steps = Int(round(tfinal / dt))
    times = [i * dt for i in 0:steps]
    cum_total_rate = 0
    threshold = -log(rand())/dt
    for _ in times[2:end]
        lambda += equation_of_motion(lambda, Nd, p)*dt
        println(lambda)
        println(Ndtotal)
        push!(lambdas, lambda)
        l = 1 + p.deltas / p.deltad * lambda
        brates = [brate_function(l, Ndi, p) for Ndi in Nd]
        urates = [urate_function(Ndi, p) for Ndi in Nd]
        total_rate = sum(brates) + sum(urates)
        cum_total_rate += total_rate
        if cum_total_rate > threshold
            event = select_event(brates, urates, total_rate)
            Nd[event[1]] += event[2]
            Ndtotal += event[2]
            cum_total_rate = 0
            threshold = -log(rand())/dt
        end
        push!(Ndtotals, Ndtotal)
    end

    return lambdas, Ndtotals, times
end
