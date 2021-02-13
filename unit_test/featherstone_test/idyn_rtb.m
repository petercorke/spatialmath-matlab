function tau = idyn_rtb(q, qd, qdd)
    twolink = [];
    mdl_twolink
    
    twolink.fast = false

    tau = twolink.rne(q, qd, qdd)
end