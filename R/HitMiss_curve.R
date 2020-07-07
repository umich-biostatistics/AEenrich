
#------------------------------------------------------------------------------------------------------------------
##------ HitMiss_Curve function: calculate miss and hit value in the KS method(considering ties) (called by get_ES function)
###
###  Input: 
###       position (vector): the index of end of each tied group
###                 eg. if the value is 11112233334, the position should be 4 6 10 11
###       hit_ind (vector): if hit, 1; Otherwise, 0
###       miss_ind (vector): if miss, 1; Otherwise, 0
###
###  Output:
###       a list with two vectors: P_hit and P_miss
#---------------------------------------------------------------------------------------------------------------
HitMiss_Curve = function(position, hit_ind, miss_ind){
  # different positions 
  n_pos = length(position)
  
  N_hit = sum(hit_ind)
  N_miss = sum(miss_ind)
  if (N_hit == 0){
    hit_value = rep(0,n_pos)
  }else{
    hit_value = cumsum(hit_ind/N_hit)
    hit_value = hit_value[position]
    #hit_value = cumsum(sapply(position, function(x) sum(which(hit_ind) %in% x)/N_hit))
  }
  miss_value = cumsum(miss_ind/N_miss)
  miss_value = miss_value[position]
  #miss_value = cumsum(sapply(position, function(x) sum(which(miss_ind) %in% x)/N_miss))
  return(list(hit=hit_value, miss=miss_value))
}
