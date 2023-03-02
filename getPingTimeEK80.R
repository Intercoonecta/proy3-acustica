# Función para extraer el 'ping time' del netCDF
# convertido desde archivos raw (EK80).
# Requiere el nombre del archivo nc con todo y ruta.
#
# 02/03/2023    
# v0.1 
#
# Ejemplo:
# pt <- getPingTimeEK80(nc = "./converted/ebes202210-D20221011-T012542.nc")
# #verificar que estan las fracciones de segundo
# strftime(pt, "%Y-%m-%d %H:%M:%OS4")

getPingTimeEK80 <- function(nc){
  require(ncdf4)
  ncf <- nc_open(nc)
  modSonar <- ncvar_get(ncf, "Sonar/sonar_software_name")
  if(modSonar[1]=="EK80") {
    ti <- ncvar_get(ncf, "Sonar/Beam_group1/ping_time")
    tu <- ncf$var$`Sonar/Beam_group1/backscatter_r`$dim[[3]]$units
    tu <-  unlist(strsplit(tu, " "))
    units <- tu[1]
    # origen
    or <- tu[3]
    or <- paste(substr(or, 1, 10), substr(or, 12, 19))
    pingTime <- as.POSIXct(ti, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS", origin = or)
    pingTime  
  } else {
    warning("Equipo no es un EK80!")
  }
  
}







