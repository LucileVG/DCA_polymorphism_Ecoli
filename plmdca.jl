using PlmDCA
using NPZ

path_msa = ARGS[1]
path_dca = ARGS[2]

println("DCA of "*path_msa)

dca_res = plmdca(path_msa, theta = 0.2, verbose = false)

J = dca_res.Jtensor                                                                                                                                                                       
h = dca_res.htensor                                                                                                                                                                       
fn_dca = dca_res.score  

npzwrite(path_dca, Dict("h_dca" => h, "J_dca" => J))

