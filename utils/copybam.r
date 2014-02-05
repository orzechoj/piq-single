bamname=rev(strsplit(bamfile,'/')[[1]])[1]
if(!file.exists(file.path(tmpdir,bamname))){
  file.copy(bamfile,file.path(tmpdir,bamname))
}
bamfile = file.path(tmpdir,bamname)