load_all_gen = function(todo=NULL){
  gen_files = list.files("../refpanel/",pattern="gen.txt.gz",full.names = T)
  legend_files = list.files("../refpanel/",pattern="legend.txt.gz", full.names = T)
  all_gen_list = list()
  for(i in 1:length(gen_files)){
    gen_file= fread(paste("gzcat",gen_files[i]),header=F)
    legend_file = fread(paste("gzcat",legend_files[i]),header=F)
    all_gen_list[[i]] = list(gen_file=gen_file, legend_file = legend_file)
    break
  }
  return(all_gen_list)
}