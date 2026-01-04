def oblivious_compaction_atomic(input_array):
    n = len(input_array)
    num_stages = n.bit_length() - 1 
    
    packets = []
    rank_counter = 0
    
    # --- oblivious Ranking ---
    for val in input_array:
        this_rank = rank_counter * val;
        packets.append({"v": val, "rank": this_rank})
        rank_counter = rank_counter + val

    # --- Oblivious Routing ---
    for stage in range(num_stages):
        stride = 1 << stage
        
        for block_start in range(0, n, 2 * stride):
            for k in range(stride):
                idx_top = block_start + k
                idx_bot = block_start + k + stride
                
                p_top = packets[idx_top]
                p_bot = packets[idx_bot]
                
                bit_top = (p_top["rank"] >> stage) & 1
                bit_bot = (p_bot["rank"] >> stage) & 1
                
                force_down_top = p_top["v"] & bit_top
                force_up_bot = p_bot["v"] & (bit_bot ^ 1)
                
                swap_cond = force_down_top | force_up_bot
                
                packets[idx_top], packets[idx_bot] = Oswap(packets[idx_top], packets[idx_bot], swap_cond)

    return [p["v"] for p in packets]